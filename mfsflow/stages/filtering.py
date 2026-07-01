"""
Filtering stage: FASTQ splitting, quality filtering, and barcode correction.

This module handles the first stage of the pipeline, which includes
splitting input FASTQ files into chunks, running quality filters,
performing barcode detection and correction, and preparing BAM chunks
for downstream mapping.
"""

import glob
import gzip
import math
import os
import subprocess

from mfsflow.bootstrap import run_barcode_discovery
from mfsflow.runtime import log_info
from mfsflow.path_layout import barcode_dir, config_dir


def estimate_avg_line_len(path, sample_lines=1000):
    """Estimate average line length from a FASTQ file by sampling.
    
    Args:
        path (str): Path to FASTQ file (optionally gzipped).
        sample_lines (int, optional): Number of lines to sample. Defaults to 1000.
        
    Returns:
        float: Average line length in bytes.
    """
    opener = gzip.open if path.endswith(".gz") else open
    total = 0
    n = 0
    with opener(path, "rb") as fh:
        for _ in range(sample_lines):
            line = fh.readline()
            if not line:
                break
            total += len(line)
            n += 1
    return (total / n) if n else 0.0


def run_filtering_stage(runtime, timer, run_stage_cmd, run_log):
    config = runtime.config
    project = runtime.project
    out_dir = runtime.out_dir
    num_threads = runtime.num_threads
    python_exec = runtime.python_exec
    samtools = runtime.tools.samtools
    pigz = runtime.tools.pigz
    seqkit = runtime.tools.seqkit
    exec_env = runtime.exec_env
    analysis_dir = runtime.analysis_dir
    tmp_merge_path = runtime.tmp_merge_path
    yaml_file = runtime.yaml_file
    resolve_script = runtime.resolve_script

    from mfsflow import pipeline_modules

    log_info("Starting Filtering Stage")

    f1_str = config.get("sequence_files", {}).get("file1", {}).get("name", "")
    f2_str = config.get("sequence_files", {}).get("file2", {}).get("name", "")

    fq1_files = [f.strip() for f in f1_str.split(",")] if f1_str else []
    fq2_files = [f.strip() for f in f2_str.split(",")] if f2_str else []

    if not fq1_files:
        raise ValueError("No file1 found in YAML configuration.")

    total_size_bytes = 0
    first_fq = fq1_files[0]

    for path in fq1_files:
        if path.endswith(".gz"):
            total_size_bytes += os.path.getsize(path) * 3
        else:
            total_size_bytes += os.path.getsize(path)

    avg_line_len = estimate_avg_line_len(first_fq, sample_lines=1000)
    if avg_line_len <= 0:
        raise ValueError(f"Failed to estimate average line length for {first_fq}")

    total_lines_est = total_size_bytes / avg_line_len
    planned_filter_jobs = max(1, max(1, num_threads) // 3)
    split_parts = max(1, min(num_threads, planned_filter_jobs * 2))

    lines_per_chunk = int(math.ceil(total_lines_est / split_parts))
    rem = lines_per_chunk % 4
    if rem != 0:
        lines_per_chunk += 4 - rem
    if lines_per_chunk < 4000:
        lines_per_chunk = 4000

    log_info(f"Total input estimation: {int(total_lines_est / 4)} reads.")
    log_info(f"Split config: {split_parts} chunk parts, {lines_per_chunk} lines per chunk.")

    with timer.section("Filtering: split FASTQ", f"parts={split_parts};lines_per_chunk={lines_per_chunk};threads={num_threads}"):
        chunk_suffixes = pipeline_modules.split_fastq(
            fq1_files,
            num_threads,
            lines_per_chunk,
            tmp_merge_path,
            project,
            pigz,
            seqkit,
            fq2_files,
            False,
            split_parts,
        )

    log_info("Running fqfilter.py on chunks")
    max_reads = config.get("counting_opts", {}).get("max_reads", 0)

    with timer.section("Filtering: fqfilter chunks", f"chunks={len(chunk_suffixes)}"):
        processes = []
        max_filter_jobs = max(1, min(len(chunk_suffixes), max(1, num_threads // 3)))
        threads_per_filter = max(1, num_threads // max_filter_jobs)
        fqfilter_pigz_threads = max(1, min(2, threads_per_filter))
        log_info(
            "fqfilter parallel jobs: "
            f"{max_filter_jobs}; pigz threads/job: {fqfilter_pigz_threads}"
        )

        def wait_for_filter_process(proc):
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"fqfilter failed (rc={proc.returncode}). Check {runtime.log_path} for details.")

        for suffix in chunk_suffixes:
            cmd = [python_exec, resolve_script("fqfilter.py"), yaml_file, pigz, suffix]
            cmd.extend(["--pigz-threads", str(fqfilter_pigz_threads)])

            if max_reads and int(max_reads) > 0:
                chunk_limit = int(int(max_reads) / len(chunk_suffixes))
                if chunk_limit < 1:
                    chunk_limit = 1
                cmd.extend(["--limit", str(chunk_limit)])

            processes.append(subprocess.Popen(cmd, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))
            if len(processes) >= max_filter_jobs:
                wait_for_filter_process(processes.pop(0))

        for proc in processes:
            wait_for_filter_process(proc)

    log_info("Cleaning up temporary FASTQ chunks...")
    with timer.section("Filtering: cleanup FASTQ chunks"):
        cleanup_candidates = (
            glob.glob(os.path.join(tmp_merge_path, "*.part_*"))
            + glob.glob(os.path.join(tmp_merge_path, "*.part_*.gz"))
            + glob.glob(os.path.join(tmp_merge_path, "*.fq.part_*"))
            + glob.glob(os.path.join(tmp_merge_path, "*.fq.part_*.gz"))
            + glob.glob(os.path.join(tmp_merge_path, "*.fastq.part_*"))
            + glob.glob(os.path.join(tmp_merge_path, "*.fastq.part_*.gz"))
        )
        removed_chunks = 0
        for path in cleanup_candidates:
            if not os.path.exists(path):
                continue
            base = os.path.basename(path)
            if base.endswith(".bam") or base.endswith(".bai") or base.endswith(".txt"):
                continue
            if ".raw.tagged." in base or ".filtered.tagged." in base:
                continue
            os.remove(path)
            removed_chunks += 1
        log_info(f"Removed temporary FASTQ chunks: {removed_chunks}")

    log_info("Merging BAM Stats")
    with timer.section("Filtering: merge BAM stats"):
        pipeline_modules.merge_bam_stats(tmp_merge_path, project, analysis_dir, yaml_file, samtools)

    bc_bin_table = os.path.join(barcode_dir(analysis_dir), f"{project}.BCbinning.txt")
    expect_id_barcode_file = os.path.join(config_dir(out_dir), "expect_id_barcode.tsv")

    if config.get("barcode_source") == "samplesheet_barcode":
        log_info("Samplesheet barcode mode: skipping barcode detection/binning")
        bc_bin_for_correction = os.devnull
    else:
        if config.get("sample", {}).get("sample_type", "").lower() == "discover":
            log_info("Running Barcode Discovery")
            with timer.section("Filtering: barcode discovery"):
                run_barcode_discovery(config, project, analysis_dir)
        log_info("Running Barcode Detection")
        run_stage_cmd([python_exec, resolve_script("barcode_detection.py"), yaml_file], "BCdetection")
        bc_bin_for_correction = bc_bin_table

    umi_chunks = []
    int_chunks = []
    if config.get("barcode_source") == "samplesheet_barcode" or os.path.exists(bc_bin_table):
        stream_bc_correction = bool(config.get("performance_opts", {}).get("stream_bc_correction", True))
        if stream_bc_correction:
            log_info("Using streaming BC correction during Mapping")
            with timer.section("Filtering: prepare raw BAM chunks for streaming correction", f"chunks={len(chunk_suffixes)}"):
                raw_chunks = []
                for suffix in chunk_suffixes:
                    raw_bam = os.path.join(tmp_merge_path, f"{project}{suffix}.raw.tagged.bam")
                    if not os.path.exists(raw_bam):
                        raise FileNotFoundError(f"Expected raw tagged chunk not found: {raw_bam}")
                    raw_chunks.append(raw_bam)
                umi_chunks = list(raw_chunks)
                int_chunks = list(raw_chunks)
        else:
            log_info("Correcting BC Tags")
            with timer.section("Filtering: correct BC tags", f"chunks={len(chunk_suffixes)}"):
                correct_processes = []

                for suffix in chunk_suffixes:
                    raw_bam = os.path.join(tmp_merge_path, f"{project}{suffix}.raw.tagged.bam")
                    fixed_bam_umi = os.path.join(tmp_merge_path, f"{project}{suffix}.filtered.tagged.umi.bam")
                    fixed_bam_int = os.path.join(tmp_merge_path, f"{project}{suffix}.filtered.tagged.internal.bam")

                    umi_chunks.append(fixed_bam_umi)
                    int_chunks.append(fixed_bam_int)

                    if os.path.exists(fixed_bam_umi):
                        os.remove(fixed_bam_umi)
                    if os.path.exists(fixed_bam_int):
                        os.remove(fixed_bam_int)

                    cmd_args = [
                        python_exec,
                        resolve_script("correct_BCtag.py"),
                        raw_bam,
                        fixed_bam_umi,
                        fixed_bam_int,
                        bc_bin_for_correction,
                        expect_id_barcode_file,
                    ]
                    correct_processes.append(subprocess.Popen(cmd_args, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))

                for proc in correct_processes:
                    proc.wait()
                    if proc.returncode != 0:
                        raise RuntimeError(f"correct_BCtag failed (rc={proc.returncode}). Check {runtime.log_path} for details.")

                for suffix in chunk_suffixes:
                    raw_bam = os.path.join(tmp_merge_path, f"{project}{suffix}.raw.tagged.bam")
                    if os.path.exists(raw_bam):
                        os.remove(raw_bam)

        log_info("Skipping physical merge of chunks (will stream to STAR)...")

    return umi_chunks, int_chunks
