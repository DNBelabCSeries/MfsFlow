"""
Counting stage: featureCounts gene quantification and DGE analysis.

This module handles the third stage of the pipeline, which performs
gene-level quantification using featureCounts and generates digital
gene expression (DGE) matrices from aligned BAM files.
"""

import os

from mfsflow.runtime import remove_path, log_info


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

    # Keep the aligned BAMs until all counting work succeeds so a failed DGE
    # step can be resumed from Counting without rerunning Mapping.
    remove_path(umi_aligned)
    remove_path(int_aligned)
    remove_path(umi_to_tx)
    remove_path(int_to_tx)

    gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    stats_enabled = str(config.get("make_stats", "yes")).lower() in ["yes", "true"]
    if not stats_enabled:
        remove_path(gene_tagged_bam)
