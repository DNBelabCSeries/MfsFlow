"""
Pipeline stage orchestration: loads configuration and executes stages sequentially.

This module provides the main pipeline runner that loads the YAML
configuration, initializes the runtime environment, and executes
the pipeline stages (Filtering, Mapping, Counting, Summarising)
in the correct order based on the user-specified starting stage.
"""

import os
import sys

import yaml

from mfsflow.commands import run_stage_cmd as run_timed_stage_cmd
from mfsflow.config import persist_run_config
from mfsflow.logging_utils import Tee, log_info
from mfsflow.preflight import run_preflight
from mfsflow.runtime import PipelineRuntime
from mfsflow.stage_state import invalidate_stage_success, record_stage_success
from mfsflow.stages import COUNTING, FILTERING, MAPPING, STAGE_ORDER, SUMMARISING
from mfsflow.stages.counting import cleanup_counting_inputs, run_counting_stage
from mfsflow.stages.filtering import run_filtering_stage
from mfsflow.stages.mapping import cleanup_mapping_inputs, run_mapping_stage
from mfsflow.stages.statistics import cleanup_statistics_inputs, run_statistics_stage
from mfsflow.timer import PipelineTimer


def run_pipeline_stages(yaml_file):
    """Orchestrate pipeline stages from a generated run_config.yaml.
    
    Loads the configuration, initializes the runtime environment, and
    executes pipeline stages sequentially starting from the specified stage.
    
    Args:
        yaml_file (str): Path to the run configuration YAML file.
    """
    log_info(f"Loading config from {yaml_file}...")
    with open(yaml_file, "r") as handle:
        config = yaml.safe_load(handle)

    runtime = PipelineRuntime.from_config(config, yaml_file)
    which_stage = runtime.which_stage
    exec_env = runtime.exec_env
    log_path = runtime.log_path
    timing_path = runtime.timing_path

    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    original_stdout = sys.stdout
    with open(log_path, "a") as run_log:
        sys.stdout = Tee(original_stdout, run_log)
        try:
            timer = PipelineTimer(timing_path, runtime.project)
            log_info("Running preflight checks...")
            with timer.section("Preflight"):
                run_preflight(config, runtime)
            persist_run_config(runtime.config, runtime.yaml_file)
            log_info("Preflight checks passed.")

            start_index = STAGE_ORDER.index(which_stage)
            for stage in STAGE_ORDER[start_index:]:
                invalidate_stage_success(runtime, stage)

            log_info(f"Starting Pipeline for project: {runtime.project}")
            log_info(f"Stage: {which_stage}")
            log_info(f"Timing log: {timing_path}")
            log_info(f"Temporary chunk directory: {runtime.tmp_merge_path}")

            def run_stage_cmd(cmd, stage_name, shell=False):
                run_timed_stage_cmd(cmd, stage_name, run_log, exec_env, timer, log_path, shell=shell)

            umi_chunks = None
            int_chunks = None

            if which_stage == "Filtering":
                umi_chunks, int_chunks = run_filtering_stage(runtime, timer, run_stage_cmd, run_log)
                persist_run_config(runtime.config, runtime.yaml_file)
                manifest = record_stage_success(runtime, FILTERING)
                log_info(f"Filtering artifact manifest: {manifest}")

            if which_stage in ["Filtering", "Mapping"]:
                mapping_umi_chunks, mapping_int_chunks = run_mapping_stage(
                    runtime,
                    run_stage_cmd,
                    umi_chunks=umi_chunks,
                    int_chunks=int_chunks,
                )
                manifest = record_stage_success(runtime, MAPPING)
                log_info(f"Mapping artifact manifest: {manifest}")
                cleanup_mapping_inputs(mapping_umi_chunks, mapping_int_chunks)

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                run_counting_stage(runtime, run_stage_cmd)
                manifest = record_stage_success(runtime, COUNTING)
                log_info(f"Counting artifact manifest: {manifest}")

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                run_statistics_stage(runtime, run_stage_cmd)
                manifest = record_stage_success(runtime, SUMMARISING)
                log_info(f"Summarising artifact manifest: {manifest}")
                # Counting inputs (Mapping BAMs) are only removed once the full
                # pipeline has reached Summarising, so Counting can still be
                # resumed before the run is fully complete.
                cleanup_counting_inputs(runtime)
                cleanup_statistics_inputs(runtime)

            log_info("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout
