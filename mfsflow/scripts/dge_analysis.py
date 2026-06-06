#!/usr/bin/env python3
"""
Digital gene expression (DGE) analysis for standalone script execution.

This module performs digital gene expression analysis, including count matrix
generation, UMI deduplication, and expression quantification for single-cell
RNA sequencing data processing.
"""

import sys
import os
import yaml
import subprocess
import csv
import gzip
import shutil
from collections import defaultdict, Counter
import scipy.io
import scipy.sparse
import numpy as np
import pandas as pd
import multiprocessing
import json
import shutil

try:
    import pysam
except ImportError:
    print("Error: 'pysam' module is required. Please install it via pip install pysam")
    sys.exit(1)

try:
    from mfsflow.scripts.umi_utils import cluster_umis
    from mfsflow.scripts.h5ad_export import export_h5ad
    from mfsflow.scripts.path_layout import barcode_dir, expression_dir, stats_dir
except ImportError:
    from umi_utils import cluster_umis
    from h5ad_export import export_h5ad
    from path_layout import barcode_dir, expression_dir, stats_dir

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def process_barcode_worker(args):
    """
    Worker function for parallel UMI clustering.
    Args:
        args: tuple (bc, umis_exon_map, umis_intron_map, ham_dist)
    Returns:
        bc, res_umi_counts_exon, res_umi_counts_intron, res_umi_counts_inex, res_correction
    """
    bc, umis_exon_map, umis_intron_map, ham_dist, need_correction_map, collect_histograms = args
    
    res_umi_counts_exon = {}
    res_umi_counts_intron = {}
    res_umi_counts_inex = {}
    res_correction = {}
    
    genes_exon = set(umis_exon_map.keys())
    genes_intron = set(umis_intron_map.keys())
    all_genes = genes_exon | genes_intron
    
    dist_counts = []
    for gene in all_genes:
        umis_ex_counts = umis_exon_map.get(gene)
        umis_in_counts = umis_intron_map.get(gene)

        if not umis_ex_counts and not umis_in_counts:
            continue

        umis_total_counts = Counter()
        if umis_ex_counts:
            umis_total_counts.update(umis_ex_counts)
        if umis_in_counts:
            umis_total_counts.update(umis_in_counts)

        mapping = cluster_umis(umis_total_counts, threshold=ham_dist)

        if ham_dist > 0 and need_correction_map:
            res_correction[gene] = mapping

        unique_total = set(mapping[u] for u in umis_total_counts.keys())
        res_umi_counts_inex[gene] = len(unique_total)

        if collect_histograms:
            canonical_counts = Counter()
            for u, cnt in umis_total_counts.items():
                canonical_counts[mapping[u]] += cnt
            dist_counts.extend(list(canonical_counts.values()))

        if umis_ex_counts:
            unique_ex = set(mapping[u] for u in umis_ex_counts.keys())
            res_umi_counts_exon[gene] = len(unique_ex)

        if umis_in_counts:
            unique_in = set(mapping[u] for u in umis_in_counts.keys())
            res_umi_counts_intron[gene] = len(unique_in)

    return bc, res_umi_counts_exon, res_umi_counts_intron, res_umi_counts_inex, res_correction, dist_counts
                    
def count_worker(args):
    """
    Worker for Pass 1: Count UMIs and Reads in a genomic region.
    Optimized for memory efficiency using compact data structures.
    """
    bam_file, chroms, barcode_set, gene_set, count_introns, collect_global_umis = args
    
    local_read_counts_raw = {
        'exon': defaultdict(lambda: defaultdict(int)),
        'intron': defaultdict(lambda: defaultdict(int))
    }
    local_umi_data = {
        'exon': defaultdict(lambda: defaultdict(Counter)),
        'intron': defaultdict(lambda: defaultdict(Counter))
    }
    local_global_umi_counts = defaultdict(Counter) if collect_global_umis else None
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for chrom in chroms:
                try:
                    iter_reads = bam.fetch(chrom)
                except ValueError:
                    continue 
                    
                for read in iter_reads:
                    if read.is_unmapped:
                        continue
                    if not read.has_tag("CB"):
                        continue
                    
                    bc = read.get_tag("CB")
                    if bc not in barcode_set:
                        continue
                    
                    if collect_global_umis and read.has_tag("UR"):
                        umi = read.get_tag("UR")
                        local_global_umi_counts[bc][umi] += 1
                    
                    try:
                        gene_id = read.get_tag("GX")
                    except KeyError:
                        continue 
                    
                    if gene_set and gene_id not in gene_set:
                        continue
                    
                    ftype = "exon"
                    if read.has_tag("RE"):
                        xf = read.get_tag("RE")
                        if xf == "N":
                            ftype = "intron"
                        elif xf == "E":
                            ftype = "exon"
                        else:
                            continue
                    
                    if not count_introns and ftype == 'intron':
                        continue
                    
                    local_read_counts_raw[ftype][bc][gene_id] += 1
                    
                    if read.has_tag("UR"):
                        umi = read.get_tag("UR")
                        local_umi_data[ftype][bc][gene_id][umi] += 1
                        
    except Exception as e:
        print(f"Error in count_worker for {chroms}: {e}")
        return None

    ret_read_counts = {
        'exon': {k: dict(v) for k, v in local_read_counts_raw['exon'].items()},
        'intron': {k: dict(v) for k, v in local_read_counts_raw['intron'].items()}
    }
    ret_umi_data = {
        'exon': {k: dict(v) for k, v in local_umi_data['exon'].items()},
        'intron': {k: dict(v) for k, v in local_umi_data['intron'].items()}
    }
    ret_global_umi = {k: dict(v) for k, v in local_global_umi_counts.items()} if collect_global_umis else {}

    return ret_read_counts, ret_umi_data, ret_global_umi

def correction_worker(args):
    """
    Worker for Pass 2: Correct UMIs.
    Writes to a temporary sorted BAM file (chunks).
    Input BAM is already sorted, so we just read and write in order.
    """
    bam_file, chroms, correction_map, ham_dist, out_tmp_bam = args
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as infile:
            with pysam.AlignmentFile(out_tmp_bam, "wb", template=infile) as outfile:
                for chrom in chroms:
                    try:
                        iter_reads = infile.fetch(chrom)
                    except ValueError:
                        continue
                        
                    for read in iter_reads:
                        raw_umi = None
                        if read.has_tag("UR"):
                            raw_umi = read.get_tag("UR")
                        elif read.has_tag("UB"):
                            raw_umi = read.get_tag("UB")
                        if not raw_umi:
                            outfile.write(read)
                            continue
                        if ham_dist <= 0:
                            read.set_tag("UB", raw_umi)
                            outfile.write(read)
                            continue

                        if read.has_tag("CB") and read.has_tag("GX"):
                            bc = read.get_tag("CB")
                            gene = read.get_tag("GX")
                            final_umi = raw_umi
                            
                            if bc in correction_map:
                                bc_map = correction_map[bc]
                                if gene in bc_map:
                                    final_umi = bc_map[gene].get(raw_umi, raw_umi)
                            
                            # Ensure final_umi is a string
                            if final_umi is None: 
                                final_umi = raw_umi
                                    
                            read.set_tag("UB", str(final_umi))
                        else:
                            # Revert: If no cell/gene assignment, remove UB tag (as requested)
                            read.set_tag("UB", None)
                        
                        outfile.write(read)
        
        return out_tmp_bam
        
    except Exception as e:
        print(f"Error in correction_worker for {chroms}: {e}")
        return None


def resolve_corrected_umi(read, correction_map, ham_dist):
    raw_umi = None
    if read.has_tag("UR"):
        raw_umi = read.get_tag("UR")
    elif read.has_tag("UB"):
        raw_umi = read.get_tag("UB")
    if not raw_umi:
        return None, False

    if ham_dist <= 0:
        return raw_umi, True

    if read.has_tag("CB") and read.has_tag("GX"):
        bc = read.get_tag("CB")
        gene = read.get_tag("GX")
        bc_map = correction_map.get(bc)
        if bc_map:
            gene_map = bc_map.get(gene)
            if gene_map:
                final_umi = gene_map.get(raw_umi, raw_umi)
                return (raw_umi if final_umi is None else str(final_umi)), True
    return raw_umi, False


def write_corrected_bam_single_pass(input_bam, out_bam, correction_map, ham_dist, threads):
    write_threads = max(1, int(threads))
    with pysam.AlignmentFile(input_bam, "rb", threads=max(1, write_threads // 2)) as infile:
        with pysam.AlignmentFile(out_bam, "wb", template=infile, threads=write_threads) as outfile:
            for read in infile.fetch(until_eof=True):
                final_umi, keep_ub = resolve_corrected_umi(read, correction_map, ham_dist)
                if final_umi is None:
                    outfile.write(read)
                    continue
                if keep_ub:
                    read.set_tag("UB", final_umi)
                else:
                    read.set_tag("UB", None)
                outfile.write(read)

def natural_sort_key(s):
    """
    Key for natural sorting (e.g., A2 < A10).
    Splits string into mixed list of strings and integers.
    """
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def load_barcodes(out_dir, project):
    """
    Loads barcodes defining the matrix columns.
    Priority 1: expect_id_barcode.tsv (Well IDs) - Ensures all wells are present & consistent.
    Priority 2: kept_barcodes.txt (Detected Seqs) - Fallback for droplet/unstructured data.
    """
    config_dir = os.path.join(out_dir, "config")
    expect_file = os.path.join(config_dir, "expect_id_barcode.tsv")
    
    barcodes = []
    source = ""
    
    if os.path.exists(expect_file):
        print(f"Loading reference barcodes (Well IDs) from {expect_file}...")
        source = "expect"
        with open(expect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # Skip header if present (check commonly used header names)
                if parts[0] in ['wellID', 'WellID', 'CellID', 'Barcode']:
                    continue
                if parts[0]:
                    barcodes.append(parts[0])
    else:
        # 2. Fallback to Kept Barcodes (Analysis Dir)
        kept_file = os.path.join(barcode_dir(out_dir), f"{project}kept_barcodes.txt")
        print(f"Loading reference barcodes (Detected) from {kept_file}...")
        source = "kept"
        if os.path.exists(kept_file):
            with open(kept_file, 'r') as f:
                # Check header
                first = f.readline()
                if not ('XC' in first or 'n' in first):
                    p = first.replace(',', '\t').split('\t')
                    if p[0]: barcodes.append(p[0])

                for line in f:
                    p = line.replace(',', '\t').split('\t')
                    if p[0]: barcodes.append(p[0])
        else:
            print("Warning: No barcode file found. Matrix will be empty.")
            return [], set()

    # Remove duplicates just in case
    barcodes = list(set(barcodes))
    
    # Sort
    # Use natural sort for consistency (e.g. P1A2 before P1A10)
    barcodes.sort(key=natural_sort_key)
    
    print(f"Loaded {len(barcodes)} barcodes from {source}.")
    return barcodes, set(barcodes)

def load_genes_from_gtf(gtf_file):
    print(f"Loading reference genes from {gtf_file}...")
    gene_order = []
    gene_map = {}
    seen_genes = set()
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue # Use exons to find gene entries
            
            attributes = parts[8]
            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            elif 'gene_id' in attributes: # Fallback
                try: gene_id = attributes.split('gene_id')[1].strip().split(';')[0].strip('"')
                except: pass
            
            if not gene_id or gene_id in seen_genes: continue
            
            gene_name = gene_id
            if 'gene_name "' in attributes:
                gene_name = attributes.split('gene_name "')[1].split('"')[0]
            
            seen_genes.add(gene_id)
            gene_order.append(gene_id)
            gene_map[gene_id] = gene_name
            
    return gene_order, gene_map

def write_sparse_matrix(counts_dict, gene_list, gene_names_map, barcode_list, out_dir, subdir_name):
    full_out_dir = os.path.join(expression_dir(out_dir), subdir_name)
    print(f"Generating deterministic, sorted sparse matrix in {full_out_dir}...")
    
    if not os.path.exists(full_out_dir):
        os.makedirs(full_out_dir)
    
    # Use fixed lists for indices
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    bc_to_idx = {b: i for i, b in enumerate(barcode_list)}
    
    rows = []
    cols = []
    data = []

    for bc, gene_counts in counts_dict.items():
        bc_idx = bc_to_idx.get(bc)
        if bc_idx is None:
            continue
        for gene, count in gene_counts.items():
            if count <= 0:
                continue
            gene_idx = gene_to_idx.get(gene)
            if gene_idx is None:
                continue
            cols.append(bc_idx)
            rows.append(gene_idx)
            data.append(int(count))

    if data:
        rows = np.asarray(rows, dtype=np.int32)
        cols = np.asarray(cols, dtype=np.int32)
        data = np.asarray(data, dtype=np.int32)
        order = np.lexsort((rows, cols))
        rows = rows[order]
        cols = cols[order]
        data = data[order]
        nnz = int(data.size)
    else:
        rows = np.asarray([], dtype=np.int32)
        cols = np.asarray([], dtype=np.int32)
        data = np.asarray([], dtype=np.int32)
        nnz = 0

    matrix_file_gz = os.path.join(full_out_dir, "matrix.mtx.gz")
    
    # Manually write MatrixMarket format to gzip stream
    # Header: %%MatrixMarket matrix coordinate integer general
    # Size: Rows Cols Entries
    with gzip.open(matrix_file_gz, 'wt') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("%\n")
        f.write(f"{len(gene_list)} {len(barcode_list)} {nnz}\n")
        
        chunk_size = 100_000
        for start in range(0, nnz, chunk_size):
            end = min(start + chunk_size, nnz)
            lines = [f"{int(rows[i]) + 1} {int(cols[i]) + 1} {int(data[i])}\n" for i in range(start, end)]
            f.write("".join(lines))

    with gzip.open(os.path.join(full_out_dir, "barcodes.tsv.gz"), "wt") as f:
        for bc in barcode_list:
            f.write(f"{bc}\n")
            
    with gzip.open(os.path.join(full_out_dir, "features.tsv.gz"), "wt") as f:
        for g_id in gene_list:
            g_name = gene_names_map.get(g_id, g_id)
            f.write(f"{g_id}\t{g_name}\tGene Expression\n")

def cluster_with_global(bc_args):
    """
    Worker function for parallel UMI clustering (includes global saturation stats).
    Moved to top-level to allow pickling.
    """
    bc, ex_map, in_map, global_counts, h_dist, need_correction_map, collect_global_umis, collect_histograms = bc_args
    # Cluster for genes
    res_bc, res_ex, res_in, res_inex, res_corr, gene_dist_counts = process_barcode_worker(
        (bc, ex_map, in_map, h_dist, need_correction_map, collect_histograms)
    )
    
    # Cluster for global saturation
    global_dist_counts = []
    if collect_global_umis and global_counts:
        mapping = cluster_umis(global_counts, threshold=h_dist)
        canonical_counts = Counter()
        for u, cnt in global_counts.items():
            canonical_counts[mapping[u]] += cnt
        global_dist_counts = list(canonical_counts.values())
    
    return res_bc, res_ex, res_in, res_inex, res_corr, global_dist_counts, gene_dist_counts

def process_bam_and_matrix(bam_file, out_bam, config, threads):
    project = config['project']
    out_dir = config['out_dir']
    ham_dist = int(config['counting_opts'].get('Ham_Dist', 0))
    count_introns = config.get('counting_opts', {}).get('introns', True)
    make_stats = bool(config.get('make_stats', True))
    make_sorted_bam = bool(config.get('make_sorted_bam', False))
    make_ub_bam = bool(config.get('make_ub_bam', False))
    need_correction_map = bool(ham_dist > 0 and (make_sorted_bam or make_ub_bam))
    collect_global_umis = bool(make_stats)
    samtools_exec = sys.argv[2] # Passed from main
    
    # Load Reference Lists
    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Updated barcode loading logic
    barcode_list, barcode_set = load_barcodes(out_dir, project)
        
    gene_list, gene_names_ref = load_genes_from_gtf(gtf_file)
    gene_set = set(gene_list)
    
    print(f"Reference: {len(barcode_list)} Barcodes (Cols), {len(gene_list)} Genes (Rows)")
    
    # Ensure BAM Index for Parallel Access
    temp_sorted_bam = None
    try:
        if not os.path.exists(bam_file + ".bai"):
            print(f"Indexing BAM {bam_file} for parallel processing...")
            pysam.index(bam_file)
    except Exception as e:
        print(f"BAM Indexing failed ({e}). Assuming BAM is not coordinate sorted.")
        print(f"Sorting input BAM to temporary file using {threads} threads...")
        temp_sorted_bam = bam_file + ".temp_sorted.bam"
        sort_cmd = [samtools_exec, "sort", "-@", str(threads), "-o", temp_sorted_bam, bam_file]
        subprocess.check_call(sort_cmd)
        print("Indexing temporary sorted BAM...")
        pysam.index(temp_sorted_bam)
        bam_file = temp_sorted_bam

    try:
        # Get Chromosomes
        with pysam.AlignmentFile(bam_file, "rb") as b:
            references = b.references
            
        # Split chromosomes into chunks for workers
        worker_count = max(1, min(int(threads), len(references) or 1))
        chunk_size = max(1, (len(references) + worker_count - 1) // worker_count)
        ref_chunks = [references[i:i + chunk_size] for i in range(0, len(references), chunk_size)]
        
        print(f"Pass 1: Parallel Counting ({worker_count} workers, {len(ref_chunks)} chunks)...")
        
        # --- PASS 1: Parallel Counting ---
        read_counts_raw = {
            'exon': defaultdict(lambda: defaultdict(int)),
            'intron': defaultdict(lambda: defaultdict(int))
        }
        umi_data = {
            'exon': defaultdict(lambda: defaultdict(Counter)),
            'intron': defaultdict(lambda: defaultdict(Counter))
        }
        global_umi_raw = defaultdict(Counter) if collect_global_umis else {}
        
        pass1_args = [
            (bam_file, chunk, barcode_set, gene_set, count_introns, collect_global_umis)
            for chunk in ref_chunks
        ]
        
        with multiprocessing.Pool(worker_count) as pool:
            for res in pool.imap_unordered(count_worker, pass1_args):
                if not res: continue
                
                partial_read, partial_umi, partial_global = res
                
                # Merge Global UMI (for saturation)
                if collect_global_umis:
                    for bc, umis in partial_global.items():
                        global_umi_raw[bc].update(umis)

                # Merge logic (In-memory reduce)
                for ftype in ['exon', 'intron']:
                    for bc, genes in partial_read[ftype].items():
                        for gene, count in genes.items():
                            read_counts_raw[ftype][bc][gene] += count
                    
                    for bc, genes in partial_umi[ftype].items():
                        for gene, umis in genes.items():
                            umi_data[ftype][bc][gene].update(umis)

        print("Pass 1 Complete. Calculating Statistics...")

        # Container for Saturation Distribution (Frequency of Read Counts)
        global_umi_freq = Counter() if collect_global_umis else None
        gene_umi_freq = Counter() if collect_global_umis else None

        # Final count containers
        final_umi_counts = {
            'exon': defaultdict(lambda: defaultdict(int)),
            'intron': defaultdict(lambda: defaultdict(int)),
            'inex': defaultdict(lambda: defaultdict(int))
        }
        
        final_read_counts = {
            'exon': read_counts_raw['exon'],
            'intron': read_counts_raw['intron'],
            'inex': defaultdict(lambda: defaultdict(int))
        }

        # Calculate Inex Read Counts
        for bc in read_counts_raw['exon']:
            for gene in read_counts_raw['exon'][bc]:
                final_read_counts['inex'][bc][gene] += read_counts_raw['exon'][bc][gene]
                
        for bc in read_counts_raw['intron']:
            for gene in read_counts_raw['intron'][bc]:
                final_read_counts['inex'][bc][gene] += read_counts_raw['intron'][bc][gene]

        # Calculate UMI Counts (with Clustering)
        all_bcs_umis = list(set(umi_data['exon'].keys()) | set(umi_data['intron'].keys()) | set(global_umi_raw.keys()))
        total_bcs = len(all_bcs_umis)
        
        print(
            f"Clustering UMIs for {total_bcs} barcodes using {worker_count} workers "
            f"(Ham_Dist={ham_dist}, correction_map={'on' if need_correction_map else 'off'}, "
            f"saturation={'on' if collect_global_umis else 'off'})..."
        )
        
        # Modified worker to also cluster global UMIs
        # cluster_with_global is now defined at module level

        pool_args = []
        for bc in all_bcs_umis:
            pool_args.append((
                bc,
                umi_data['exon'].get(bc, {}),
                umi_data['intron'].get(bc, {}),
                global_umi_raw.get(bc, {}) if collect_global_umis else {},
                ham_dist,
                need_correction_map,
                collect_global_umis,
                collect_global_umis,
            ))
            
        correction_map = defaultdict(lambda: defaultdict(dict)) if need_correction_map else {}
        
        with multiprocessing.Pool(worker_count) as pool:
             for i, res in enumerate(pool.imap_unordered(cluster_with_global, pool_args, chunksize=20)):
                if i % 100 == 0: print(f"Clustering {i}/{total_bcs}...", end='\r')
                
                bc, c_ex, c_in, c_inex, corr, g_dist, gene_dist = res
                
                if c_ex: final_umi_counts['exon'][bc].update(c_ex)
                if c_in: final_umi_counts['intron'][bc].update(c_in)
                if c_inex: final_umi_counts['inex'][bc].update(c_inex)
                if need_correction_map and corr:
                    correction_map[bc].update(corr)
                if collect_global_umis and g_dist:
                    global_umi_freq.update(g_dist)
                if collect_global_umis and gene_dist:
                    gene_umi_freq.update(gene_dist)

        print("\nWriting Matrices...")
        
        # Calculate and print correction stats
        if need_correction_map:
            total_corrections = 0
            for bc in correction_map:
                for gene in correction_map[bc]:
                    for child, parent in correction_map[bc][gene].items():
                        if child != parent:
                            total_corrections += 1
            print(f"Total UMI corrections found: {total_corrections}")
            if ham_dist > 0 and total_corrections == 0:
                print("WARNING: Ham_Dist > 0 but no corrections found. Check input data or barcode matching.")

        
        # Write Saturation Data (Histogram)
        if collect_global_umis:
            sat_file = os.path.join(stats_dir(out_dir), f"{project}.saturation_dist.json")
            gene_sat_file = os.path.join(stats_dir(out_dir), f"{project}.gene_saturation_dist.json")
            if not os.path.exists(os.path.dirname(sat_file)):
                os.makedirs(os.path.dirname(sat_file))

            print(f"Writing saturation distribution histograms to {os.path.dirname(sat_file)}...")
            with open(sat_file, 'w') as f:
                json.dump(dict(global_umi_freq), f)
            with open(gene_sat_file, 'w') as f:
                json.dump(dict(gene_umi_freq), f)

        # Write Matrices
        write_sparse_matrix(final_umi_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.umi")
        if count_introns:
            write_sparse_matrix(final_umi_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.umi")
            write_sparse_matrix(final_umi_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.umi")
            
        write_sparse_matrix(final_read_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.read")
        if count_introns:
            write_sparse_matrix(final_read_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.read")
            write_sparse_matrix(final_read_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.read")

        if bool(config.get('make_h5ad', True)):
            print("Exporting combined H5AD...")
            h5ad_path = export_h5ad(out_dir, project, config=config)
            print(f"H5AD written: {h5ad_path}")

        if not make_sorted_bam and not make_ub_bam:
            return

        if make_ub_bam and not make_sorted_bam:
            if not out_bam:
                out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.bam")
            print(f"Writing UB-corrected BAM in a single pass ({threads} threads for BAM IO)...")
            write_corrected_bam_single_pass(bam_file, out_bam, correction_map, ham_dist, threads)
            return

        # Input BAM is coordinate-sorted here either because it already had a
        # valid index or because we created a temporary sorted copy above, so a
        # single-pass write preserves sorted order without per-chrom chunk BAMs.
        print(f"Writing UB-corrected sorted BAM in a single pass ({threads} threads for BAM IO)...")
        write_corrected_bam_single_pass(bam_file, out_bam, correction_map, ham_dist, threads)

        print("Indexing Final BAM...")
        pysam.index(out_bam)
        
    finally:
        if temp_sorted_bam and os.path.exists(temp_sorted_bam):
            print(f"Cleaning up temporary sorted BAM: {temp_sorted_bam}")
            os.remove(temp_sorted_bam)
            if os.path.exists(temp_sorted_bam + ".bai"):
                os.remove(temp_sorted_bam + ".bai")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 dge_analysis.py <yaml_config> <samtools_exec>")
        sys.exit(1)
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 1))
    input_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    if not os.path.exists(input_bam):
        print(f"Error: Input BAM {input_bam} not found.")
        sys.exit(1)
    make_sorted_bam = bool(config.get('make_sorted_bam', False))
    make_ub_bam = bool(config.get('make_ub_bam', False))
    out_bam = None
    if make_sorted_bam:
        out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
    elif make_ub_bam:
        out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.bam")
    process_bam_and_matrix(input_bam, out_bam, config, threads=num_threads)
    print("DGE Analysis pipeline finished.")

if __name__ == "__main__":
    main()
