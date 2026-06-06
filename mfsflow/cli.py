"""
Command-line interface for the MfsFlow analysis pipeline.

This module provides the argument parser and main entry point for the
MfsFlow analysis pipeline, handling user input, configuration building,
and pipeline execution.
"""

import argparse
import os
import sys

from mfsflow.stages import FILTERING, STAGE_ORDER


def build_parser():
    """Build the command-line argument parser for the MfsFlow pipeline.
    
    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="MfsFlow Data Analysis Pipeline")
    parser.add_argument("--fastqs", required=True, help="Directory containing input R1/R2 FASTQ files")
    parser.add_argument("--samplesheet", help="CSV samplesheet for equal-length R1/R2 data")
    parser.add_argument("--genomeDir", required=True, help="Reference directory containing star/ and genes/genes.gtf or genes.gtf.gz")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--outdir", help="Output directory (default: ./<sample_name>)")
    parser.add_argument("--threads", type=int, default=20, help="Number of threads")
    parser.add_argument("--tmpRoot", help="Temporary root for chunk BAM/FASTQ files, e.g. /dev/shm")
    parser.add_argument("--stage", choices=STAGE_ORDER, default=FILTERING, help="Analysis stage to start from")

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--manual", help='Manual sample IDs (comma separated, e.g. "20,21"). Sets sample_type=manual.')
    mode_group.add_argument("--plate", help='Plate ID (e.g. "1"). Sets sample_type=auto.')
    mode_group.add_argument("--expectBarcode", help="Path to custom barcode file. Sets sample_type=custom.")
    return parser


def generate_report(config, analysis_failed=False):
    """Generate the HTML report for the completed pipeline analysis.
    
    Args:
        config (dict): Pipeline configuration dictionary containing project
            settings and output directory paths.
    """
    from mfsflow.path_layout import logs_dir
    from mfsflow.runtime import PipelineTimer

    from mfsflow import report
    from mfsflow.runtime import log_info

    log_info('Generating HTML Report...')
    report_config = dict(config)
    report_config["_analysis_failed"] = bool(analysis_failed)
    timing_path = os.path.join(logs_dir(config["out_dir"]), "pipeline_timing.tsv")
    report_timer = PipelineTimer(timing_path, config["project"])
    with report_timer.section("Report: HTML generation"):
        return report.generate_multi_report(config["project"], config["out_dir"], report_config)


def main(argv=None):
    """Main entry point for the MfsFlow analysis pipeline.
    
    Parses command-line arguments, builds configuration, runs pipeline stages,
    and generates the final HTML report.
    
    Args:
        argv (list, optional): Command-line arguments. Defaults to sys.argv.
    """
    args = build_parser().parse_args(argv)

    from mfsflow.preflight import check_python_dependencies

    check_python_dependencies()

    import time
    from pathlib import Path

    from mfsflow.pipeline_config import build_base_config, resolve_samplesheet_barcodes
    from mfsflow.run_config import write_run_config
    from mfsflow.bootstrap import create_barcode_tables, create_output_dirs, process_fastq_inputs
    from mfsflow.config.validation import require_supported_python, validate_input_files
    from mfsflow.pipeline.runner import run_pipeline_stages
    from mfsflow.runtime import format_duration, log_info

    require_supported_python()

    log_info(f'Start analysis for {args.sample}.')

    script_dir = str(Path(__file__).resolve().parent)
    config, samplesheet_records = build_base_config(args, script_dir)

    create_output_dirs(config)
    log_info('Directories created.')

    validate_input_files(config)

    process_fastq_inputs(config)
    log_info('Fastq processed.')

    create_barcode_tables(config)
    log_info('Barcode files created.')

    if config.get("barcode_source") == "samplesheet_barcode":
        expect_id_file = os.path.join(config["out_dir"], "config", "expect_id_barcode.tsv")
        config["fastq_groups"] = resolve_samplesheet_barcodes(samplesheet_records, expect_id_file)

    final_yaml_path = write_run_config(config)

    log_info(f'Config generated: {final_yaml_path}')
    log_info('Starting Pipeline...')

    pipeline_start = time.perf_counter()
    try:
        run_pipeline_stages(final_yaml_path)
    except Exception:
        from mfsflow.runtime import log_error

        pipeline_duration = time.perf_counter() - pipeline_start
        log_error(
            f'Analysis failed after {format_duration(pipeline_duration)}. '
            'Attempting to generate a partial report from completed outputs.'
        )
        try:
            generate_report(config, analysis_failed=True)
        except Exception as report_exc:
            log_error(f"Partial report generation also failed: {report_exc}")
        raise
    else:
        pipeline_duration = time.perf_counter() - pipeline_start
        log_info(f'All analysis finished (Duration: {format_duration(pipeline_duration)}).')
        generate_report(config)


if __name__ == "__main__":
    main()
