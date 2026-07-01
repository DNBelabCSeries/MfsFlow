"""Stage artifact manifests, success markers, and resume validation."""

import glob
import json
import os
import shutil
import subprocess
import gzip
from datetime import datetime
from pathlib import Path

from mfsflow.path_layout import expression_dir, stage_state_dir, stats_dir
from mfsflow.stages import COUNTING, FILTERING, MAPPING, SUMMARISING


def _nonempty_file(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _artifact_record(path):
    stat = os.stat(path)
    return {
        "path": os.path.abspath(path),
        "size_bytes": stat.st_size,
        "modified_at": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def _validate_gzip(path, label):
    try:
        with gzip.open(path, "rb") as handle:
            if not handle.read(1):
                raise RuntimeError(f"{label} is empty after decompression: {path}")
    except (OSError, EOFError) as exc:
        raise RuntimeError(f"{label} is corrupt or unreadable: {path}: {exc}") from exc


def _quickcheck_bams(runtime, paths, unmapped=False):
    """Validate BAM structure, allowing header-only targets for unmapped BAMs."""
    paths = list(paths)
    if not paths:
        return
    samtools = getattr(getattr(runtime, "tools", None), "samtools", None)
    samtools_available = bool(samtools) and (
        (os.path.isfile(str(samtools)) and os.access(str(samtools), os.X_OK))
        or shutil.which(str(samtools))
    )
    if samtools_available:
        command = [str(samtools), "quickcheck", "-v"]
        if unmapped:
            command.append("-u")
        command.extend(paths)
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode == 0:
            return
        if not unmapped:
            details = result.stdout.strip() or "samtools quickcheck returned a non-zero status"
            raise RuntimeError(f"BAM integrity check failed: {details}")
        # Unmapped BAMs lack @SQ headers; samtools < 1.10 (no -u support) rejects
        # them. Fall back to pysam which handles no-SQ BAMs natively.
    if unmapped:
        _pysam_validate_bams(paths)


def _pysam_validate_bams(paths):
    """Fallback validation for unmapped (no-SQ) BAMs when samtools is unavailable
    or lacks quickcheck -u support (samtools < 1.10)."""
    import pysam
    failures = []
    for path in paths:
        try:
            with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
                next(bam.fetch(until_eof=True), None)
        except Exception as exc:
            failures.append(f"{path}: {exc}")
    if failures:
        raise RuntimeError("BAM integrity check failed: " + "; ".join(failures))


def _glob_files(patterns):
    found = set()
    for pattern in patterns:
        found.update(path for path in glob.glob(pattern) if _nonempty_file(path))
    return sorted(found)


def stage_artifacts(runtime, stage):
    """Return key non-empty artifacts produced by a completed stage."""
    project = runtime.project
    out_dir = runtime.out_dir
    tmp_dir = runtime.tmp_merge_path

    if stage == FILTERING:
        return _glob_files([
            os.path.join(tmp_dir, f"{project}*.raw.tagged.bam"),
            os.path.join(tmp_dir, f"{project}*.filtered.tagged.umi.bam"),
            os.path.join(tmp_dir, f"{project}*.filtered.tagged.internal.bam"),
            os.path.join(out_dir, "barcodes", f"{project}.*"),
            os.path.join(out_dir, "stats", f"{project}.Q30stats.txt"),
            os.path.join(out_dir, "config", "expect_id_barcode.tsv"),
        ])
    if stage == MAPPING:
        return _glob_files([
            os.path.join(out_dir, f"{project}.filtered.tagged.*.Aligned.out.bam"),
            os.path.join(out_dir, f"{project}.filtered.tagged.*.Aligned.toTranscriptome.out.bam"),
        ])
    if stage == COUNTING:
        return _glob_files([
            os.path.join(expression_dir(out_dir), f"{project}.*", "matrix.mtx.gz"),
            os.path.join(expression_dir(out_dir), f"{project}.*", "barcodes.tsv.gz"),
            os.path.join(expression_dir(out_dir), f"{project}.*", "features.tsv.gz"),
            os.path.join(expression_dir(out_dir), f"{project}.h5ad"),
            os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged*.bam"),
            os.path.join(stats_dir(out_dir), f"{project}.*saturation_dist.json"),
        ])
    if stage == SUMMARISING:
        return _glob_files([
            os.path.join(stats_dir(out_dir), f"{project}.stats.tsv"),
            os.path.join(stats_dir(out_dir), f"{project}.saturation.tsv"),
            os.path.join(stats_dir(out_dir), f"{project}.geneBodyCoverage.txt"),
            os.path.join(stats_dir(out_dir), f"{project}.read_stats.json"),
        ])
    raise ValueError(f"Unknown stage: {stage}")


def _validate_stage_outputs(runtime, stage, artifacts):
    project = runtime.project
    out_dir = runtime.out_dir

    if stage == FILTERING:
        has_mapping_input = any(
            name.endswith((".raw.tagged.bam", ".filtered.tagged.umi.bam"))
            for name in artifacts
        )
        if not has_mapping_input:
            raise RuntimeError("Filtering completed but no non-empty Mapping input BAM chunks were found.")
        _quickcheck_bams(runtime, [path for path in artifacts if path.endswith(".bam")], unmapped=True)
    elif stage == MAPPING:
        umi_bam = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
        if not _nonempty_file(umi_bam):
            raise RuntimeError(f"Mapping completed but required UMI aligned BAM is missing or empty: {umi_bam}")
        _quickcheck_bams(runtime, [path for path in artifacts if path.endswith(".bam")])
    elif stage == COUNTING:
        matrix = os.path.join(expression_dir(out_dir), f"{project}.exon.umi", "matrix.mtx.gz")
        if not _nonempty_file(matrix):
            raise RuntimeError(f"Counting completed but required expression matrix is missing or empty: {matrix}")
        _validate_gzip(matrix, "Expression matrix")
    elif stage == SUMMARISING and str(runtime.config.get("make_stats", "yes")).lower() in ("yes", "true"):
        table = os.path.join(stats_dir(out_dir), f"{project}.stats.tsv")
        if not _nonempty_file(table):
            raise RuntimeError(f"Summarising completed but required stats table is missing or empty: {table}")


def invalidate_stage_success(runtime, stage):
    """Remove stale success state before executing a stage."""
    state_dir = stage_state_dir(runtime.out_dir)
    for suffix in (".success", ".manifest.json"):
        path = os.path.join(state_dir, f"{stage}{suffix}")
        if os.path.exists(path):
            os.remove(path)


def record_stage_success(runtime, stage):
    """Validate outputs and atomically write a manifest and success marker."""
    artifacts = stage_artifacts(runtime, stage)
    _validate_stage_outputs(runtime, stage, artifacts)

    state_dir = Path(stage_state_dir(runtime.out_dir))
    state_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = state_dir / f"{stage}.manifest.json"
    marker_path = state_dir / f"{stage}.success"
    manifest_tmp = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
    marker_tmp = marker_path.with_suffix(marker_path.suffix + ".tmp")
    payload = {
        "project": runtime.project,
        "stage": stage,
        "status": "success",
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "artifacts": [_artifact_record(path) for path in artifacts],
    }
    try:
        manifest_tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(manifest_tmp, manifest_path)
        marker_tmp.write_text(payload["completed_at"] + "\n", encoding="utf-8")
        os.replace(marker_tmp, marker_path)
    finally:
        for path in (manifest_tmp, marker_tmp):
            if path.exists():
                path.unlink()
    return manifest_path


def validate_resume_inputs(runtime):
    """Validate inputs required to start from the configured resume stage."""
    stage = runtime.which_stage
    project = runtime.project
    out_dir = runtime.out_dir

    if stage == FILTERING:
        return
    if stage == MAPPING:
        candidates = _glob_files([
            os.path.join(runtime.tmp_merge_path, f"{project}*.raw.tagged.bam"),
            os.path.join(runtime.tmp_merge_path, f"{project}*.filtered.tagged.umi.bam"),
            os.path.join(out_dir, f"{project}.filtered.tagged.umi.unmapped.bam"),
        ])
        if not candidates:
            raise RuntimeError(
                "Cannot resume from Mapping: no non-empty Filtering BAM artifacts were found. "
                "Resume from Filtering or restore the Filtering outputs."
            )
        _quickcheck_bams(runtime, candidates, unmapped=True)
    elif stage == COUNTING:
        required = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
        if not _nonempty_file(required):
            raise RuntimeError(
                f"Cannot resume from Counting: required Mapping artifact is missing or empty: {required}"
            )
        _quickcheck_bams(runtime, [required])
    elif stage == SUMMARISING:
        required = os.path.join(expression_dir(out_dir), f"{project}.exon.umi", "matrix.mtx.gz")
        if not _nonempty_file(required):
            raise RuntimeError(
                f"Cannot resume from Summarising: required Counting artifact is missing or empty: {required}"
            )
        _validate_gzip(required, "Counting expression matrix")
