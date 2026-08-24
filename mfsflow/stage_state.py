"""Stage artifact manifests, success markers, and resume validation."""

import glob
import json
import os
import shutil
import subprocess
import gzip
import hashlib
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
        "mtime_ns": stat.st_mtime_ns,
        "modified_at": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def _config_digest(config, digest_version=2):
    stable_config = json.loads(json.dumps(config, default=str))
    # The requested resume stage changes between runs and is not an analysis
    # parameter, so it must not invalidate an earlier stage manifest.
    stable_config.pop("which_Stage", None)
    stable_config.pop("which_stage", None)
    sample_config = stable_config.get("sample")
    if isinstance(sample_config, dict):
        for key in list(sample_config):
            if key.startswith("discovered_"):
                sample_config.pop(key, None)
    stable_config.pop("_analysis_failed", None)
    stable_config.pop("tool_versions", None)
    if digest_version >= 2:
        # Execution placement and parallelism may change during a resume without
        # changing the analytical inputs or requested outputs.
        stable_config.pop("num_threads", None)
        stable_config.pop("toolkit_directory", None)
        performance = stable_config.get("performance_opts")
        if isinstance(performance, dict):
            for key in ("tmp_root", "tool_cache", "min_free_gb", "disk_space_multiplier"):
                performance.pop(key, None)
    encoded = json.dumps(stable_config, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


# Compressed-size threshold (bytes) for the lightweight gzip check performed
# immediately after Counting. Resume validation always requests a full stream
# read, regardless of file size.
_GZIP_FULL_CHECK_MAX_BYTES = 64 * 1024 * 1024


def _validate_gzip(path, label, full_check=True):
    try:
        size = os.path.getsize(path)
        if size == 0:
            raise RuntimeError(f"{label} is empty: {path}")
        with gzip.open(path, "rb") as handle:
            # Always confirm we can actually decompress data from the stream.
            if not handle.read(1):
                raise RuntimeError(f"{label} is empty after decompression: {path}")
            if full_check or size <= _GZIP_FULL_CHECK_MAX_BYTES:
                # Small file: read the whole stream to catch mid-stream corruption.
                while handle.read(1024 * 1024):
                    pass
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
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
    _pysam_validate_bams(paths)


def _pysam_validate_bams(paths):
    """Fallback validation for unmapped (no-SQ) BAMs when samtools is unavailable
    or lacks quickcheck -u support (samtools < 1.10)."""
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "BAM integrity validation requires either executable samtools or the pysam package."
        ) from exc
    failures = []
    for path in paths:
        try:
            with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
                check_truncation = getattr(bam, "check_truncation", None)
                if check_truncation is not None:
                    check_truncation()
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
            os.path.join(out_dir, "stats", f"{project}.q30_stats.tsv"),
            os.path.join(out_dir, f"{project}.BCstats.txt"),
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
            os.path.join(stats_dir(out_dir), f"{project}.cell_matrix_stats.json"),
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
        _validate_gzip(matrix, "Expression matrix", full_check=False)
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
        "config_digest_version": 2,
        "config_sha256": _config_digest(runtime.config),
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


def validate_stage_manifest(runtime, stage):
    """Validate a previously completed stage and all recorded artifacts."""
    state_dir = stage_state_dir(runtime.out_dir)
    manifest_path = os.path.join(state_dir, f"{stage}.manifest.json")
    marker_path = os.path.join(state_dir, f"{stage}.success")
    if not os.path.isfile(manifest_path) or not os.path.isfile(marker_path):
        raise RuntimeError(f"Stage {stage} has no complete success manifest.")

    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Stage {stage} manifest is unreadable: {manifest_path}: {exc}") from exc

    if payload.get("status") != "success" or payload.get("stage") != stage:
        raise RuntimeError(f"Stage {stage} manifest is not marked successful: {manifest_path}")
    if payload.get("project") != runtime.project:
        raise RuntimeError(f"Stage {stage} manifest belongs to a different project: {manifest_path}")
    expected_digest = payload.get("config_sha256")
    digest_version = int(payload.get("config_digest_version", 1))
    if expected_digest and expected_digest != _config_digest(runtime.config, digest_version=digest_version):
        raise RuntimeError(f"Stage {stage} manifest was generated from a different configuration.")

    artifacts = payload.get("artifacts") or []
    if not artifacts:
        if stage == SUMMARISING and str(runtime.config.get("make_stats", "yes")).lower() not in ("yes", "true", "1"):
            return payload
        raise RuntimeError(f"Stage {stage} manifest contains no artifacts: {manifest_path}")
    changed = []
    for artifact in artifacts:
        path = artifact.get("path", "")
        if not _nonempty_file(path):
            changed.append(f"missing or empty: {path}")
            continue
        expected_size = artifact.get("size_bytes")
        if expected_size is not None and os.path.getsize(path) != expected_size:
            changed.append(f"size changed: {path}")
            continue
        expected_mtime = artifact.get("mtime_ns")
        if expected_mtime is not None and os.stat(path).st_mtime_ns != expected_mtime:
            changed.append(f"modified after manifest: {path}")
    if changed:
        hint = ""
        if stage in (FILTERING, MAPPING, COUNTING):
            hint = (
                " These are intermediate artifacts that are removed once the pipeline "
                "reaches the next stages. If you are resuming after a fully completed run, "
                "re-run from the Filtering stage (--stage Filtering) to regenerate them."
            )
        raise RuntimeError(
            f"Stage {stage} artifacts failed manifest validation: " + "; ".join(changed) + hint
        )

    if stage == FILTERING:
        _quickcheck_bams(runtime, [a["path"] for a in artifacts if a.get("path", "").endswith(".bam")], unmapped=True)
    elif stage == MAPPING:
        _quickcheck_bams(runtime, [a["path"] for a in artifacts if a.get("path", "").endswith(".bam")])
    elif stage == COUNTING:
        for artifact in artifacts:
            path = artifact.get("path", "")
            if path.endswith(".mtx.gz"):
                _validate_gzip(path, "Expression matrix", full_check=True)
    return payload


def validate_resume_inputs(runtime):
    """Validate inputs required to start from the configured resume stage."""
    stage = runtime.which_stage
    project = runtime.project
    out_dir = runtime.out_dir

    if stage == FILTERING:
        return
    if stage == MAPPING:
        if os.path.exists(os.path.join(stage_state_dir(out_dir), f"{FILTERING}.manifest.json")):
            validate_stage_manifest(runtime, FILTERING)
            return
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
        if os.path.exists(os.path.join(stage_state_dir(out_dir), f"{MAPPING}.manifest.json")):
            validate_stage_manifest(runtime, MAPPING)
            return
        required = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
        if not _nonempty_file(required):
            raise RuntimeError(
                f"Cannot resume from Counting: required Mapping artifact is missing or empty: {required}"
            )
        _quickcheck_bams(runtime, [required])
    elif stage == SUMMARISING:
        if os.path.exists(os.path.join(stage_state_dir(out_dir), f"{COUNTING}.manifest.json")):
            validate_stage_manifest(runtime, COUNTING)
            return
        required = os.path.join(expression_dir(out_dir), f"{project}.exon.umi", "matrix.mtx.gz")
        if not _nonempty_file(required):
            raise RuntimeError(
                f"Cannot resume from Summarising: required Counting artifact is missing or empty: {required}"
            )
        _validate_gzip(required, "Counting expression matrix")
