"""
Counting stage: featureCounts gene quantification and DGE analysis.

This module handles the third stage of the pipeline, which performs
gene-level quantification using featureCounts and generates digital
gene expression (DGE) matrices from aligned BAM files.
"""

import glob
import os
import shutil

from mfsflow.fs_utils import remove_path
from mfsflow.logging_utils import log_info
from mfsflow.path_layout import expression_dir, stats_dir


def _clear_previous_counting_outputs(runtime):
    """Remove Counting and downstream outputs while preserving Mapping BAMs."""
    project = runtime.project
    analysis_dir = runtime.analysis_dir
    candidates = []
    candidates.extend(glob.glob(os.path.join(expression_dir(analysis_dir), f"{project}.*")))
    candidates.extend(glob.glob(os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged*")))
    candidates.extend(
        path
        for path in glob.glob(os.path.join(analysis_dir, f"{project}.filtered.tagged.*.Aligned.out.bam.*"))
        if not path.endswith((".bai", ".csi"))
    )

    stats_root = stats_dir(analysis_dir)
    for suffix in (
        ".read_stats.json",
        ".saturation_dist.json",
        ".gene_saturation_dist.json",
        ".cell_matrix_stats.json",
        ".stats.tsv",
        ".saturation.tsv",
        ".geneBodyCoverage.txt",
        ".geneBodyCoverage.pdf",
        ".features.pdf",
    ):
        candidates.append(os.path.join(stats_root, project + suffix))

    removed = 0
    for path in dict.fromkeys(candidates):
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
            removed += 1
        elif os.path.isfile(path) or os.path.islink(path):
            os.unlink(path)
            removed += 1
    if removed:
        log_info(f"Removed {removed} stale Counting/Summarising artifact(s) before rerun.")


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
    _clear_previous_counting_outputs(runtime)

    umi_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
    int_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")

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
