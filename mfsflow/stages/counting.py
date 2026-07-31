"""
Counting stage: featureCounts gene quantification and DGE analysis.

This module handles the third stage of the pipeline, which performs
gene-level quantification using featureCounts and generates digital
gene expression (DGE) matrices from aligned BAM files.
"""

import os

from mfsflow.fs_utils import remove_path
from mfsflow.logging_utils import log_info


def run_counting_stage(runtime, run_stage_cmd):
    """Execute the counting stage of the pipeline.
    
    Args:
        runtime (PipelineRuntime): Pipeline runtime configuration.
        run_stage_cmd (callable): Function to run stage commands with timing.
    """
    project = runtime.project
    analysis_dir = runtime.analysis_dir
    yaml_file = runtime.yaml_file
    python_exec = runtime.python_exec
    samtools = runtime.tools.samtools
    resolve_script = runtime.resolve_script
    config = runtime.config

    log_info("Starting Counting Stage")

    umi_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
    int_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
    umi_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.toTranscriptome.out.bam")
    int_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.toTranscriptome.out.bam")

    featurecounts_cmd = [
        python_exec,
        resolve_script("run_featurecounts.py"),
        yaml_file,
        "--umi_bam",
        umi_aligned,
        "--internal_bam",
        int_aligned,
    ]
    run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

    log_info("Starting DGE Analysis (Python)")
    dge_cmd = [python_exec, resolve_script("dge_analysis.py"), yaml_file, samtools]
    run_stage_cmd(dge_cmd, "dge_analysis.py")

    gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    stats_enabled = str(config.get("make_stats", "yes")).lower() in ["yes", "true"]
    if not stats_enabled:
        remove_path(gene_tagged_bam)


def cleanup_counting_inputs(runtime):
    """Remove Mapping BAMs only after Counting has been marked successful."""
    project = runtime.project
    for suffix in (
        ".filtered.tagged.umi.Aligned.out.bam",
        ".filtered.tagged.internal.Aligned.out.bam",
        ".filtered.tagged.umi.Aligned.toTranscriptome.out.bam",
        ".filtered.tagged.internal.Aligned.toTranscriptome.out.bam",
    ):
        path = os.path.join(runtime.analysis_dir, project + suffix)
        try:
            remove_path(path)
        except OSError as exc:
            log_info(f"Warning: could not remove Counting input {path}: {exc}")
