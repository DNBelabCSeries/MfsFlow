#!/usr/bin/env python3
"""
Mapping analysis for standalone script execution.

This module performs read mapping analysis, including alignment statistics,
mapping quality assessment, and read distribution analysis for single-cell
RNA sequencing data processing.
"""

import sys
import os
import glob
import subprocess
import yaml
import math
import shutil
import collections
import itertools
import gzip
import sys

try:
    from mfsflow.scripts.path_layout import barcode_dir
except ImportError:
    from path_layout import barcode_dir

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def run_cmd(cmd, shell=True):
    print(f"Running: {cmd}")
    subprocess.check_call(cmd, shell=shell)

def get_bam_read_length(bam_file, samtools, n_reads=1000):
    proc = subprocess.Popen(
        [samtools, 'view', bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    lengths = collections.Counter()
    try:
        for line in itertools.islice(proc.stdout, n_reads):
            parts = line.rstrip('\n').split('\t')
            if len(parts) > 9:
                lengths[len(parts[9])] += 1
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        proc.terminate()
        proc.wait()

    # rc=0: success
    # rc=-15: SIGTERM (we terminated it)
    # rc=-13: SIGPIPE (we closed the pipe while it was writing, expected)
    if proc.returncode not in (0, -15, -13):
        stderr = proc.stderr.read().strip() if proc.stderr else ''
        raise RuntimeError(f"samtools view failed (rc={proc.returncode}) for {bam_file}: {stderr}")

    if not lengths:
        return 0

    return lengths.most_common(1)[0][0]


def get_stream_corrected_read_length(bam_files, bc_bin_file, expect_id_file, target_type, n_reads=1000):
    """
    Estimate read length after the same streaming barcode correction used before STAR.
    This preserves STAR sjdbOverhang behavior when raw tagged chunks are streamed
    directly instead of pre-writing corrected UMI/internal BAM chunks.
    """
    try:
        import pysam
        try:
            from mfsflow.scripts.barcode_corrector import correct_read_barcode, load_bc_map, load_id_map
        except ImportError:
            from barcode_corrector import correct_read_barcode, load_bc_map, load_id_map
    except ImportError:
        return 0

    bc_map = load_bc_map(bc_bin_file)
    id_map, internal_bcs = load_id_map(expect_id_file, strict=True)
    lengths = collections.Counter()
    seen = 0

    for bam_file in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
        except Exception:
            continue
        with bam:
            for read in bam:
                correction = correct_read_barcode(read, bc_map, id_map, internal_bcs)
                if correction is None:
                    if target_type != "umi":
                        continue
                else:
                    if target_type == "umi" and correction.is_internal:
                        continue
                    if target_type == "internal" and not correction.is_internal:
                        continue

                query_sequence = read.query_sequence
                if query_sequence:
                    lengths[len(query_sequence)] += 1
                    seen += 1
                    if seen >= n_reads:
                        return lengths.most_common(1)[0][0]

    if not lengths:
        return 0
    return lengths.most_common(1)[0][0]


def determine_target_read_length(bam_files, target_type, streams_raw_chunks, bc_bin_file, expect_id_file, samtools):
    read_len = 0
    if streams_raw_chunks and bam_files:
        read_len = get_stream_corrected_read_length(
            bam_files,
            bc_bin_file,
            expect_id_file,
            target_type,
        )
    if read_len <= 0 and bam_files:
        read_len = get_bam_read_length(bam_files[0], samtools)
    return read_len


def get_dir_size_gb(path):
    try:
        out = subprocess.check_output(['du', '-sk', path], text=True).split('\t', 1)[0]
        kb = float(out.strip())
        return kb / (1024 * 1024)
    except Exception:
        return 25.0


def read_star_index_params(star_index):
    params = {}
    param_file = os.path.join(star_index, "genomeParameters.txt")
    if not os.path.exists(param_file):
        return params

    with open(param_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if "\t" in line:
                key, value = line.split("\t", 1)
            else:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                key, value = parts
            params[key.strip()] = value.strip()
    return params


def read_star_index_chromosomes(star_index):
    chr_file = os.path.join(star_index, "chrName.txt")
    if not os.path.exists(chr_file):
        return set()
    with open(chr_file, "r") as handle:
        return {line.strip() for line in handle if line.strip()}


def star_index_has_embedded_sjdb(index_params):
    for key in ("sjdbGTFfile", "sjdbFileChrStartEnd", "sjdbInsertSave"):
        value = index_params.get(key, "")
        if value and value not in ("-", "None", "0"):
            return True
    overhang = index_params.get("sjdbOverhang", "")
    return bool(overhang and overhang not in ("-", "None", "0"))


def build_star_misc_base(star_index, num_threads, read_layout, index_has_sjdb, final_gtf, read_len):
    misc_parts = [
        f"--genomeDir {star_index}",
        f"--runThreadN {num_threads}",
        f"--readFilesType SAM {read_layout}",
        "--outSAMmultNmax 1",
        "--outFilterMultimapNmax 50",
        "--outSAMunmapped Within",
        "--outSAMtype BAM Unsorted",
        "--limitOutSJcollapsed 5000000",
    ]
    if not index_has_sjdb:
        if read_len <= 0:
            raise ValueError("Unable to determine read length required for STAR sjdbOverhang.")
        misc_parts.insert(1, f"--sjdbGTFfile {final_gtf}")
        misc_parts.insert(2, f"--sjdbOverhang {read_len - 1}")
    return " ".join(misc_parts)

def setup_gtf(config, project, out_dir, samtools):
    # Handle additional files
    gtf = config['reference']['GTF_file']
    additional_files = config['reference'].get('additional_files', [])
    if isinstance(additional_files, str):
        additional_files = [additional_files] if additional_files.strip() else []
    star_index = config['reference'].get('STAR_index', '')
    if not additional_files:
        final_gtf = os.path.join(out_dir, f"{project}.final_annot.gtf")
        if str(gtf).endswith(".gz"):
            with gzip.open(gtf, "rt") as infile, open(final_gtf, "w") as outfile:
                shutil.copyfileobj(infile, outfile)
        else:
            shutil.copyfile(gtf, final_gtf)
        return final_gtf, ""
    
    # Process additional fasta
    add_gtf_path = os.path.join(out_dir, "additional_sequence_annot.gtf")
    additional_chroms = set()
    with open(add_gtf_path, 'w') as out_f:
        for fa in additional_files:
            # Get lengths using samtools faidx
            # Assuming faidx exists? If not generate it.
            if not os.path.exists(fa + ".fai"):
                subprocess.run([samtools, "faidx", fa], check=True)
            
            with open(fa + ".fai") as fai:
                for line in fai:
                    parts = line.split('\t')
                    name = parts[0]
                    length = parts[1]
                    additional_chroms.add(name)
                    # Write custom GTF line
                    # R code: gene_id "name"; transcript_id "name"; ...
                    attr = f'gene_id "{name}"; transcript_id "{name}"; exon_number "1"; gene_name "{name}"; gene_biotype "User"; transcript_name "{name}"; exon_id "{name}"'
                    out_f.write(f"{name}\tUser\texon\t1\t{length}\t.\t+\t.\t{attr}\n")

    index_chroms = read_star_index_chromosomes(star_index)
    missing_from_index = sorted(additional_chroms - index_chroms)
    if missing_from_index:
        raise ValueError(
            "additional_files contigs must already be included in the STAR index. "
            "Rebuild STAR_index before mapping. Missing contigs: "
            + ", ".join(missing_from_index)
        )

    final_gtf = os.path.join(out_dir, f"{project}.final_annot.gtf")
    with open(final_gtf, 'w') as outfile:
        # Cat original GTF
        opener = gzip.open if str(gtf).endswith(".gz") else open
        with opener(gtf, 'rt') as infile:
            shutil.copyfileobj(infile, outfile)
        # Cat additional
        with open(add_gtf_path, 'r') as infile:
            shutil.copyfileobj(infile, outfile)
            
    return final_gtf, ""

def run_star_pipe(corrector_args, star_cmd_str):
    """
    Runs a pipeline: Corrector (Python) -> STAR
    """
    print(f"Starting Pipeline: {' '.join(corrector_args[:3])}... -> STAR")
    
    # Start Producer (Corrector)
    # Use list args to avoid shell quoting issues with many files
    p1 = subprocess.Popen(corrector_args, stdout=subprocess.PIPE)
    
    # Start Consumer (STAR)
    # star_cmd_str should use --readFilesIn /dev/stdin
    p2 = subprocess.Popen(star_cmd_str, shell=True, stdin=p1.stdout)
    
    # Close p1's stdout in this parent process so only p2 holds it
    p1.stdout.close()
    
    # Wait for completion
    p2.wait()
    p1.wait()
    
    if p2.returncode != 0:
        raise RuntimeError(f"STAR failed with return code {p2.returncode}")
    
    # If STAR succeeds, p1 should also succeed (0) or SIGPIPE (141/-13)
    if p1.returncode not in (0, -13, 141):
         print(f"Warning: Corrector script exited with code {p1.returncode}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file')
    parser.add_argument('--umi_bam', required=False, help="Merged UMI Unmapped BAM")
    parser.add_argument('--internal_bam', required=False, help="Merged Internal Unmapped BAM")
    parser.add_argument('--expect_id_file', required=False, help="Path to expect_id_barcode.tsv")
    args = parser.parse_args()
        
    yaml_file = args.yaml_file
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 1))
    
    samtools = config.get('samtools_exec', 'samtools')
    star_exec = config.get('STAR_exec', 'STAR')
    star_index = config['reference']['STAR_index']
    
    # Executables Check
    if not shutil.which(star_exec) and not os.path.exists(star_exec):
         print(f"Error: STAR executable not found: {star_exec}")
         sys.exit(1)

    # 1. Parse Inputs
    # Support both command line args (comma separated list) and legacy/yaml lookup
    umi_bams = []
    if args.umi_bam:
        umi_bams = [x.strip() for x in args.umi_bam.split(',') if x.strip()]
    
    internal_bams = []
    if args.internal_bam:
        internal_bams = [x.strip() for x in args.internal_bam.split(',') if x.strip()]
        
    # Explicit expect_id_file argument overrides everything
    if args.expect_id_file:
        expect_id_file = args.expect_id_file
    else:
        # Fallback (legacy logic)
        barcode_config_path = config['barcodes'].get('barcode_file')
        if barcode_config_path:
            config_dir = os.path.dirname(barcode_config_path)
            expect_id_file = os.path.join(config_dir, "expect_id_barcode.tsv")
        else:
            root_dir = os.path.dirname(out_dir.rstrip(os.sep))
            expect_id_file = os.path.join(root_dir, "config", "expect_id_barcode.tsv")

    if not os.path.exists(expect_id_file):
        raise FileNotFoundError(f"ID Map file not found: {expect_id_file}")
    
    # Check existence
    for f in umi_bams:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input UMI BAM not found: {f}")
    for f in internal_bams:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input Internal BAM not found: {f}")

    if not umi_bams and not internal_bams:
        raise FileNotFoundError("No input unmapped BAMs found.")

    # 2. Setup GTF
    final_gtf, param_add_fa = setup_gtf(config, project, out_dir, samtools)
    star_index_params = read_star_index_params(star_index)
    index_has_sjdb = star_index_has_embedded_sjdb(star_index_params)

    # Define paths for corrector
    corrector_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stream_corrector.py")
    bc_bin_file = os.path.join(barcode_dir(out_dir), f"{project}.BCbinning.txt")
    if not os.path.exists(bc_bin_file):
        bc_bin_file = os.devnull

    streams_raw_chunks = any(os.path.basename(x).endswith(".raw.tagged.bam") for x in (umi_bams + internal_bams))
    umi_read_len = determine_target_read_length(
        umi_bams,
        "umi",
        streams_raw_chunks,
        bc_bin_file,
        expect_id_file,
        samtools,
    ) if umi_bams else 0
    internal_read_len = determine_target_read_length(
        internal_bams,
        "internal",
        streams_raw_chunks,
        bc_bin_file,
        expect_id_file,
        samtools,
    ) if internal_bams else 0

    if umi_read_len > 0:
        print(f"Detected UMI Read Length: {umi_read_len}")
    if internal_read_len > 0:
        print(f"Detected Internal Read Length: {internal_read_len}")
    if index_has_sjdb:
        print("Detected STAR index already contains sjdb/GTF annotation; skipping on-the-fly --sjdbGTFfile injection.")
    else:
        print("Detected STAR index without embedded sjdb annotation; enabling on-the-fly --sjdbGTFfile injection.")

    # 4. Resource Allocation
    print(f"Allocating {num_threads} threads for sequential execution.")

    # 5. Build STAR Commands
    read_layout = config.get('read_layout', 'SE')
    
    extra_params = config['reference'].get('additional_STAR_params', '') or ""
    
    # Two-pass mode
    twopass = ""
    if config['counting_opts'].get('twoPass', False):
        twopass = "--twopassMode Basic"

    # Run UMI
    if umi_bams:
        prefix_umi = os.path.join(out_dir, f"{project}.filtered.tagged.umi.")
        
        # Corrector Args (Producer)
        corrector_args = [sys.executable or 'python3', corrector_script, '--binning', bc_bin_file, '--idmap', expect_id_file, '--type', 'umi'] + umi_bams
        
        # STAR Command (Consumer)
        # --readFilesIn /dev/stdin
        misc_base = build_star_misc_base(star_index, num_threads, read_layout, index_has_sjdb, final_gtf, umi_read_len)
        cmd_umi = f"{star_exec} {misc_base} {extra_params} {param_add_fa} {twopass} --readFilesIn /dev/stdin --outFileNamePrefix {prefix_umi}"
        
        run_star_pipe(corrector_args, cmd_umi)
        print("STAR UMI finished.")

    # Run Internal
    if internal_bams:
        prefix_int = os.path.join(out_dir, f"{project}.filtered.tagged.internal.")
        
        # Corrector Args
        corrector_args = [sys.executable or 'python3', corrector_script, '--binning', bc_bin_file, '--idmap', expect_id_file, '--type', 'internal'] + internal_bams
        
        # STAR Command
        misc_base = build_star_misc_base(star_index, num_threads, read_layout, index_has_sjdb, final_gtf, internal_read_len)
        cmd_int = f"{star_exec} {misc_base} {extra_params} {param_add_fa} {twopass} --readFilesIn /dev/stdin --outFileNamePrefix {prefix_int}"
        
        run_star_pipe(corrector_args, cmd_int)
        print("STAR Internal finished.")
        


if __name__ == "__main__":
    main()
