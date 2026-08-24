"""
Mapping stage: STAR alignment of filtered and barcode-corrected BAM files.

This module handles the second stage of the pipeline, which aligns
UMI-tagged and internal-barcode-tagged BAM files to the reference
genome using STAR, producing aligned BAM files for downstream counting.
"""

import glob
import json
import os
import shutil

from mfsflow.fs_utils import remove_path
from mfsflow.logging_utils import log_info
from mfsflow.path_layout import stage_state_dir
from mfsflow.stages import FILTERING


def _clear_previous_mapping_outputs(runtime):
    """Remove STAR outputs without deleting Filtering inputs used for resume."""
    removed = 0
    for read_type in ("umi", "internal"):
        prefix = os.path.join(
            runtime.analysis_dir,
            f"{runtime.project}.filtered.tagged.{read_type}.",
        )
        legacy_input = prefix + "unmapped.bam"
        for path in glob.glob(prefix + "*"):
            if path in {legacy_input, legacy_input + ".bai"}:
                continue
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.unlink(path)
            removed += 1
    if removed:
        log_info(f"Removed {removed} stale Mapping artifact(s) before rerun.")


def _filtering_manifest_bams(runtime):
    """Return validated Filtering BAM paths, including an earlier tmp root."""
    manifest_path = os.path.join(stage_state_dir(runtime.out_dir), f"{FILTERING}.manifest.json")
    if not os.path.isfile(manifest_path):
        return []
    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, ValueError):
        return []
    return sorted(
        artifact.get("path")
        for artifact in (payload.get("artifacts") or [])
        if str(artifact.get("path") or "").endswith(".bam")
        and os.path.isfile(artifact.get("path") or "")
    )


def run_mapping_stage(runtime, run_stage_cmd, umi_chunks=None, int_chunks=None):
    """Execute the mapping stage of the pipeline.
    
    Args:
        runtime (PipelineRuntime): Pipeline runtime configuration.
        run_stage_cmd (callable): Function to run stage commands with timing.
        umi_chunks (list, optional): UMI BAM chunk paths.
        int_chunks (list, optional): Internal barcode BAM chunk paths.
    """
    project = runtime.project
    analysis_dir = runtime.analysis_dir
    tmp_merge_path = runtime.tmp_merge_path
    out_dir = runtime.out_dir
    python_exec = runtime.python_exec
    yaml_file = runtime.yaml_file
    resolve_script = runtime.resolve_script

    log_info("Starting Mapping Stage")
    _clear_previous_mapping_outputs(runtime)

    umi_chunks = list(umi_chunks or [])
    int_chunks = list(int_chunks or [])
    umi_arg = ""
    int_arg = ""

    def find_chunks(suffix_pattern):
        found = glob.glob(os.path.join(tmp_merge_path, suffix_pattern))
        return sorted(found)

    manifest_bams = _filtering_manifest_bams(runtime) if runtime.which_stage == "Mapping" else []

    def find_manifest_chunks(suffix):
        return [path for path in manifest_bams if path.endswith(suffix)]

    if umi_chunks:
        umi_arg = ",".join(umi_chunks)
    else:
        disk_umi_chunks = find_chunks(f"{project}*.filtered.tagged.umi.bam") or find_manifest_chunks(
            ".filtered.tagged.umi.bam"
        )
        if disk_umi_chunks:
            log_info(f"Found {len(disk_umi_chunks)} UMI chunks on disk.")
            umi_chunks = disk_umi_chunks
            umi_arg = ",".join(disk_umi_chunks)
        else:
            disk_raw_chunks = find_chunks(f"{project}*.raw.tagged.bam") or find_manifest_chunks(".raw.tagged.bam")
            if disk_raw_chunks:
                log_info(f"Found {len(disk_raw_chunks)} raw tagged chunks on disk for streaming UMI correction.")
                umi_chunks = disk_raw_chunks
                umi_arg = ",".join(disk_raw_chunks)
            else:
                legacy_umi = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                if os.path.exists(legacy_umi):
                    umi_arg = legacy_umi
                elif runtime.which_stage == "Mapping":
                    raise FileNotFoundError(
                        f"Could not find input BAMs for Mapping stage. Checked for chunks in {tmp_merge_path} and merged file {legacy_umi}"
                    )

    if int_chunks:
        int_arg = ",".join(int_chunks)
    else:
        disk_int_chunks = find_chunks(f"{project}*.filtered.tagged.internal.bam") or find_manifest_chunks(
            ".filtered.tagged.internal.bam"
        )
        if disk_int_chunks:
            log_info(f"Found {len(disk_int_chunks)} Internal chunks on disk.")
            int_chunks = disk_int_chunks
            int_arg = ",".join(disk_int_chunks)
        else:
            disk_raw_chunks = find_chunks(f"{project}*.raw.tagged.bam") or find_manifest_chunks(".raw.tagged.bam")
            if disk_raw_chunks:
                log_info(f"Found {len(disk_raw_chunks)} raw tagged chunks on disk for streaming internal correction.")
                int_chunks = disk_raw_chunks
                int_arg = ",".join(disk_raw_chunks)
            else:
                legacy_int = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                if os.path.exists(legacy_int):
                    int_arg = legacy_int

    map_cmd = [
        python_exec,
        resolve_script("mapping_analysis.py"),
        yaml_file,
        "--umi_bam",
        umi_arg,
        "--internal_bam",
        int_arg,
    ]
    expect_id_file = os.path.join(out_dir, "config", "expect_id_barcode.tsv")
    map_cmd.extend(["--expect_id_file", expect_id_file])
    run_stage_cmd(map_cmd, "mapping_analysis.py")
    return umi_chunks, int_chunks

def cleanup_mapping_inputs(umi_chunks=None, int_chunks=None):
    """Remove Filtering inputs only after Mapping has been marked successful."""
    for path in set((umi_chunks or []) + (int_chunks or [])):
        try:
            remove_path(path)
        except OSError as exc:
            log_info(f"Warning: could not remove Mapping input {path}: {exc}")
