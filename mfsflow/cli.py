"""
Command-line interface for the MfsFlow analysis pipeline.

This module provides the argument parser and main entry point for the
MfsFlow analysis pipeline, handling user input, configuration building,
and pipeline execution.
"""

import argparse
import os

from mfsflow import __version__
from mfsflow.run_lock import OutputRunLock
from mfsflow.stages import FILTERING, STAGE_ORDER


def build_parser():
    """Build the command-line argument parser for the MfsFlow pipeline.
    
    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="MfsFlow Data Analysis Pipeline")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--fastqs", help="Directory containing input R1/R2 FASTQ files")
    parser.add_argument("--samplesheet", help="CSV samplesheet for equal-length R1/R2 data")
    parser.add_argument("--genomeDir", help="Reference directory containing star/ and genes/genes.gtf or genes.gtf.gz")
    parser.add_argument("--sample", help="Sample name")
    parser.add_argument("--outdir", help="Output directory (default: ./<sample_name>)")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 20; resume keeps the original value)")
    parser.add_argument("--tmpRoot", help="Temporary root for chunk BAM/FASTQ files, e.g. /dev/shm")
    parser.add_argument("--stage", choices=STAGE_ORDER, help="Analysis stage to start from")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an existing --outdir using its saved run_config.yaml without rebuilding inputs or barcode tables",
    )

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
    from mfsflow.timer import PipelineTimer

    from mfsflow import report
    from mfsflow.logging_utils import log_info

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
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be >= 1")

    if args.resume:
        if not args.outdir:
            parser.error("--resume requires --outdir")
        if not args.stage:
            parser.error("--resume requires an explicit --stage")
        conflicting = [
            option
            for option in ("fastqs", "genomeDir", "samplesheet", "manual", "plate", "expectBarcode")
            if getattr(args, option, None)
        ]
        if conflicting:
            parser.error(
                "--resume uses saved analytical inputs; remove: "
                + ", ".join(f"--{option}" for option in conflicting)
            )
        return _resume_existing_run(args, parser)

    if args.stage and args.stage != FILTERING:
        parser.error(f"--stage {args.stage} requires --resume and an existing --outdir")

    missing = [option for option in ("fastqs", "genomeDir", "sample") if not getattr(args, option)]
    if missing:
        parser.error("the following arguments are required for a new run: " + ", ".join(f"--{name}" for name in missing))
    if args.threads is None:
        args.threads = 20
    if args.stage is None:
        args.stage = FILTERING

    from pathlib import Path

    from mfsflow.bootstrap import create_output_dirs
    from mfsflow.config import build_base_config, require_supported_python
    from mfsflow.logging_utils import log_info

    require_supported_python()

    log_info(f'Start analysis for {args.sample}.')

    script_dir = str(Path(__file__).resolve().parent)
    config, samplesheet_records = build_base_config(args, script_dir)

    create_output_dirs(config)
    log_info('Directories created.')

    with OutputRunLock(config["out_dir"]):
        _validate_output_ownership(config)
        _run_analysis(config, samplesheet_records)


def _validate_output_ownership(config):
    """Reject accidental reuse of an output directory owned by another project."""
    from mfsflow.config import load_run_config

    config_path = os.path.join(config["out_dir"], "config", "run_config.yaml")
    if not os.path.isfile(config_path):
        return
    existing = load_run_config(config_path)
    existing_project = str(existing.get("project") or "")
    requested_project = str(config.get("project") or "")
    if existing_project and existing_project != requested_project:
        raise RuntimeError(
            f"Output directory already belongs to project {existing_project!r}; "
            f"refusing to write project {requested_project!r}: {config['out_dir']}"
        )


def _resume_config_path(outdir):
    """Resolve run_config.yaml from either a project root or XPRESS_PROCESSING."""
    outdir = os.path.abspath(outdir)
    processing_dir = outdir if os.path.basename(outdir) == "XPRESS_PROCESSING" else os.path.join(outdir, "XPRESS_PROCESSING")
    return os.path.join(processing_dir, "config", "run_config.yaml")


def _resume_existing_run(args, parser):
    from mfsflow.bootstrap import create_output_dirs
    from mfsflow.config import load_run_config, persist_run_config, require_supported_python
    from mfsflow.logging_utils import log_info

    require_supported_python()
    config_path = _resume_config_path(args.outdir)
    if not os.path.isfile(config_path):
        parser.error(f"resume configuration not found: {config_path}")
    config = load_run_config(config_path)
    processing_dir = os.path.dirname(os.path.dirname(config_path))
    configured_outdir = os.path.abspath(config.get("out_dir", ""))
    if configured_outdir != processing_dir:
        parser.error(
            f"saved out_dir does not match --outdir: {configured_outdir or '<missing>'} != {processing_dir}"
        )
    if args.sample and args.sample != config.get("project"):
        parser.error(
            f"--sample {args.sample!r} does not match saved project {config.get('project')!r}"
        )

    config["which_Stage"] = args.stage
    config.pop("which_stage", None)
    if args.threads is not None:
        config["num_threads"] = args.threads
    if args.tmpRoot:
        config.setdefault("performance_opts", {})["tmp_root"] = os.path.abspath(args.tmpRoot)

    create_output_dirs(config)
    with OutputRunLock(config["out_dir"]):
        persist_run_config(config, config_path)
        log_info(f"Resuming {config['project']} from {args.stage} using {config_path}.")
        _execute_pipeline(config, config_path)


def _run_analysis(config, samplesheet_records):
    """Run setup, pipeline stages, and report generation under the output lock."""
    from mfsflow.bootstrap import create_barcode_tables
    from mfsflow.config import resolve_samplesheet_fastq_groups, validate_input_files, write_run_config
    from mfsflow.logging_utils import log_info
    from mfsflow.preflight import check_python_dependencies

    validate_input_files(config)
    check_python_dependencies(config=config)

    create_barcode_tables(config)
    log_info('Barcode files created.')

    if config.get("barcode_source") == "samplesheet_barcode":
        expect_id_file = os.path.join(config["out_dir"], "config", "expect_id_barcode.tsv")
        config["fastq_groups"] = resolve_samplesheet_fastq_groups(
            samplesheet_records,
            expect_id_file,
            config.get("sample", {}).get("sample_type"),
        )

    final_yaml_path = write_run_config(config)

    log_info(f'Config generated: {final_yaml_path}')
    log_info('Starting Pipeline...')

    _execute_pipeline(config, final_yaml_path)


def _execute_pipeline(config, final_yaml_path):
    """Execute a new or resumed pipeline and report from its persisted state."""
    import time

    from mfsflow.config import load_run_config
    from mfsflow.logging_utils import log_error, log_info
    from mfsflow.pipeline.runner import run_pipeline_stages
    from mfsflow.timer import format_duration

    def current_config():
        try:
            return load_run_config(final_yaml_path)
        except Exception:
            return config

    pipeline_start = time.perf_counter()
    try:
        run_pipeline_stages(final_yaml_path)
    except Exception as pipeline_exc:
        pipeline_duration = time.perf_counter() - pipeline_start
        log_error(
            f'Analysis failed after {format_duration(pipeline_duration)}. '
            'Attempting to generate a partial report from completed outputs.'
        )
        try:
            failed_config = current_config()
            failed_config["_analysis_error"] = str(pipeline_exc)
            generate_report(failed_config, analysis_failed=True)
        except Exception as report_exc:
            log_error(f"Partial report generation also failed: {report_exc}")
        raise
    else:
        pipeline_duration = time.perf_counter() - pipeline_start
        log_info(f'All analysis finished (Duration: {format_duration(pipeline_duration)}).')
        generate_report(current_config())


if __name__ == "__main__":
    main()
