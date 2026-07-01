"""
Statistics stage: generation of quality metrics and summary statistics.

This module handles the final stage of the pipeline, which computes
and aggregates quality statistics, generating reports and summary
metrics for the completed analysis.
"""

import os

from mfsflow.fs_utils import remove_path
from mfsflow.logging_utils import log_info


def run_statistics_stage(runtime, run_stage_cmd):
    """Execute the statistics stage of the pipeline.
    
    Args:
        runtime (PipelineRuntime): Pipeline runtime configuration.
        run_stage_cmd (callable): Function to run stage commands with timing.
    """
    config = runtime.config
    if str(config.get("make_stats", "yes")).lower() not in ["yes", "true"]:
        return

    log_info("Starting Statistics Stage")
    stats_cmd = [runtime.python_exec, runtime.resolve_script("generate_stats.py"), runtime.yaml_file]
    run_stage_cmd(stats_cmd, "Stats (Python)")

    gene_tagged_bam = os.path.join(runtime.analysis_dir, f"{runtime.project}.filtered.Aligned.GeneTagged.bam")
    remove_path(gene_tagged_bam)
