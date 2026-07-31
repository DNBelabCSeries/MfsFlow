"""
Mapping stage: STAR alignment of filtered and barcode-corrected BAM files.

This module handles the second stage of the pipeline, which aligns
UMI-tagged and internal-barcode-tagged BAM files to the reference
genome using STAR, producing aligned BAM files for downstream counting.
"""

import glob
import os

from mfsflow.fs_utils import remove_path
from mfsflow.logging_utils import log_info


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

    umi_chunks = list(umi_chunks or [])
    int_chunks = list(int_chunks or [])
    umi_arg = ""
    int_arg = ""

    def find_chunks(suffix_pattern):
        found = glob.glob(os.path.join(tmp_merge_path, suffix_pattern))
        return sorted(found)

    if umi_chunks:
        umi_arg = ",".join(umi_chunks)
    else:
        disk_umi_chunks = find_chunks(f"{project}*.filtered.tagged.umi.bam")
        if disk_umi_chunks:
            log_info(f"Found {len(disk_umi_chunks)} UMI chunks on disk.")
            umi_chunks = disk_umi_chunks
            umi_arg = ",".join(disk_umi_chunks)
        else:
            disk_raw_chunks = find_chunks(f"{project}*.raw.tagged.bam")
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
        disk_int_chunks = find_chunks(f"{project}*.filtered.tagged.internal.bam")
        if disk_int_chunks:
            log_info(f"Found {len(disk_int_chunks)} Internal chunks on disk.")
            int_chunks = disk_int_chunks
            int_arg = ",".join(disk_int_chunks)
        else:
            disk_raw_chunks = find_chunks(f"{project}*.raw.tagged.bam")
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
