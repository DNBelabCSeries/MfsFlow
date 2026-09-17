#!/usr/bin/env python3
"""
Digital gene expression (DGE) analysis for standalone script execution.

This module performs digital gene expression analysis, including count matrix
generation, UMI deduplication, and expression quantification for single-cell
RNA sequencing data processing.
"""

import sys
import os
import subprocess
import gzip
import tempfile
from collections import defaultdict, Counter
import multiprocessing
import json
import time
from concurrent.futures import ProcessPoolExecutor

try:
    import pysam
except ImportError:
    pysam = None

try:
    from mfsflow.scripts.umi_utils import cluster_umis
    from mfsflow.scripts.h5ad_export import export_h5ad
    from mfsflow.scripts.dge_utils import (
        balance_reference_chunks,
        close_pass1_store,
        bounded_results,
        finalize_pass1_store,
        load_pass1_read_bundle,
        load_pass1_umi_bundle,
        open_pass1_store,
        pass1_barcode_workloads,
        pass1_barcodes,
        resolve_worker_count,
        summarize_exon_intron_counts,
        store_pass1_result,
        workload_order,
    )
    from mfsflow.scripts.read_utils import is_pair_representative
    from mfsflow.path_layout import barcode_dir, expression_dir, stats_dir, load_config
except ImportError:
    from umi_utils import cluster_umis
    from h5ad_export import export_h5ad
    from dge_utils import (
        balance_reference_chunks,
        close_pass1_store,
        bounded_results,
        finalize_pass1_store,
        load_pass1_read_bundle,
        load_pass1_umi_bundle,
        open_pass1_store,
        pass1_barcode_workloads,
        pass1_barcodes,
        resolve_worker_count,
        summarize_exon_intron_counts,
        store_pass1_result,
        workload_order,
    )
    from read_utils import is_pair_representative
    from path_layout import barcode_dir, expression_dir, stats_dir, load_config


def _as_bool(value, default=False):
    """Interpret YAML booleans and quoted boolean values consistently."""
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(value)


def _write_json_atomic(path, payload):
    """Write a JSON artifact through a sibling temporary file."""
    directory = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, temporary_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.",
        suffix=".tmp",
        dir=directory,
    )
    os.close(fd)
    try:
        with open(temporary_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)


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
            # Keep a gene marker, but retain only actual substitutions. Reads
            # whose UMI is absent from this sparse map still need the same
            # unchanged-UB behavior in resolve_corrected_umi().
            res_correction[gene] = {
                child: parent
                for child, parent in mapping.items()
                if child != parent
            }

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
                iter_reads = bam.fetch(chrom)
                for read in iter_reads:
                    # BAM retains both PE mates, but PE read/UMI/saturation
                    # statistics use the primary R1, matching zUMIs.
                    if not is_pair_representative(read):
                        continue
                    if read.is_unmapped:
                        continue
                    # Single CB fetch: replaces has_tag + get_tag (2 calls -> 1)
                    try:
                        bc = read.get_tag("CB")
                    except KeyError:
                        continue
                    if bc not in barcode_set:
                        continue

                    # Single UR fetch: replaces up to 4 has_tag/get_tag calls below
                    try:
                        umi = read.get_tag("UR")
                    except KeyError:
                        umi = None

                    if collect_global_umis and umi is not None:
                        local_global_umi_counts[bc][umi] += 1

                    try:
                        gene_id = read.get_tag("GX")
                    except KeyError:
                        continue

                    if gene_set and gene_id not in gene_set:
                        continue

                    ftype = "exon"
                    try:
                        xf = read.get_tag("RE")
                        if xf == "N":
                            ftype = "intron"
                        elif xf != "E":
                            continue
                    except KeyError:
                        pass  # no RE tag: default exon

                    if not count_introns and ftype == 'intron':
                        continue

                    local_read_counts_raw[ftype][bc][gene_id] += 1

                    if umi is not None:
                        local_umi_data[ftype][bc][gene_id][umi] += 1
                        
    except Exception as exc:
        labels = ", ".join(chroms) or "<no references>"
        raise RuntimeError(f"DGE counting failed for references [{labels}]: {exc}") from exc

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


def resolve_corrected_umi(read, correction_map, ham_dist):
    try:
        raw_umi = read.get_tag("UR")
    except KeyError:
        try:
            raw_umi = read.get_tag("UB")
        except KeyError:
            raw_umi = None
    if not raw_umi:
        return None, False

    if ham_dist <= 0:
        return raw_umi, True

    try:
        bc = read.get_tag("CB")
        gene = read.get_tag("GX")
    except KeyError:
        return raw_umi, False

    bc_map = correction_map.get(bc)
    if bc_map is not None:
        gene_map = bc_map.get(gene)
        if gene_map is not None:
            final_umi = gene_map.get(raw_umi, raw_umi)
            return (raw_umi if final_umi is None else str(final_umi)), True
    return raw_umi, False


def write_corrected_bam_single_pass(input_bam, out_bam, correction_map, ham_dist, threads):
    write_threads = max(1, int(threads))
    output_dir = os.path.dirname(os.path.abspath(out_bam)) or "."
    os.makedirs(output_dir, exist_ok=True)
    fd, temporary_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(out_bam)}.",
        suffix=".tmp.bam",
        dir=output_dir,
    )
    os.close(fd)
    try:
        with pysam.AlignmentFile(input_bam, "rb", threads=max(1, write_threads // 2)) as infile:
            with pysam.AlignmentFile(temporary_path, "wb", template=infile, threads=write_threads) as outfile:
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
        os.replace(temporary_path, out_bam)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)

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
    
    # Use fixed indices and stream entries in the same column/row order as the
    # previous lexsort implementation. The old implementation retained three
    # Python lists plus three NumPy arrays for the entire matrix, which was a
    # significant peak for large projects.
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}

    def count_entries():
        for barcode in barcode_list:
            gene_counts = counts_dict.get(barcode)
            if not gene_counts:
                continue
            for gene, count in gene_counts.items():
                if count > 0 and gene_to_idx.get(gene) is not None:
                    yield 1

    def iter_entries():
        for col_idx, barcode in enumerate(barcode_list):
            gene_counts = counts_dict.get(barcode)
            if not gene_counts:
                continue
            entries = []
            for gene, count in gene_counts.items():
                if count <= 0:
                    continue
                row_idx = gene_to_idx.get(gene)
                if row_idx is not None:
                    entries.append((row_idx, int(count)))
            for row_idx, count in sorted(entries):
                yield row_idx + 1, col_idx + 1, count

    # Matrix Market requires nnz before the data lines. Iterate twice instead
    # of materialising all entries; the count dictionaries are already in RAM.
    # The first pass only counts entries; sorting is needed only for the
    # writing pass and should not be repeated unnecessarily.
    nnz = sum(count_entries())

    final_paths = {
        "matrix": os.path.join(full_out_dir, "matrix.mtx.gz"),
        "barcodes": os.path.join(full_out_dir, "barcodes.tsv.gz"),
        "features": os.path.join(full_out_dir, "features.tsv.gz"),
    }
    temporary_paths = {}

    def make_temp_path(final_path):
        fd, temp_path = tempfile.mkstemp(
            prefix=f".{os.path.basename(final_path)}.",
            suffix=".tmp",
            dir=full_out_dir,
        )
        os.close(fd)
        return temp_path

    try:
        # Fill this mapping incrementally so a failure while creating a later
        # temporary file still lets the finally block remove earlier files.
        for key, final_path in final_paths.items():
            temporary_paths[key] = make_temp_path(final_path)

        # Write all three files completely before publishing any of them. This
        # prevents an interrupted run from leaving a valid-looking partial MEX.
        with gzip.open(temporary_paths["matrix"], "wt") as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write("%\n")
            f.write(f"{len(gene_list)} {len(barcode_list)} {nnz}\n")

            lines = []
            for row_idx, col_idx, count in iter_entries():
                lines.append(f"{row_idx} {col_idx} {count}\n")
                if len(lines) >= 100_000:
                    f.write("".join(lines))
                    lines.clear()
            if lines:
                f.write("".join(lines))

        with gzip.open(temporary_paths["barcodes"], "wt") as f:
            for bc in barcode_list:
                f.write(f"{bc}\n")

        with gzip.open(temporary_paths["features"], "wt") as f:
            for g_id in gene_list:
                g_name = gene_names_map.get(g_id, g_id)
                f.write(f"{g_id}\t{g_name}\tGene Expression\n")

        for key, final_path in final_paths.items():
            os.replace(temporary_paths[key], final_path)
            temporary_paths[key] = None
    finally:
        for temp_path in temporary_paths.values():
            if temp_path and os.path.exists(temp_path):
                os.unlink(temp_path)

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

def process_bam_and_matrix(bam_file, out_bam, config, threads, samtools_exec=None):
    analysis_started = time.perf_counter()

    def log_phase(label, started):
        print(f"DGE timing: {label} ({time.perf_counter() - started:.2f}s)", flush=True)

    project = config['project']
    out_dir = config['out_dir']
    ham_dist = int(config['counting_opts'].get('Ham_Dist', 0))
    count_introns = _as_bool(config.get('counting_opts', {}).get('introns', True), default=True)
    make_stats = _as_bool(config.get('make_stats', True), default=True)
    make_sorted_bam = _as_bool(config.get('make_sorted_bam', False))
    make_ub_bam = _as_bool(config.get('make_ub_bam', False))
    need_correction_map = bool(ham_dist > 0 and (make_sorted_bam or make_ub_bam))
    collect_global_umis = bool(make_stats)
    if samtools_exec is None:
        samtools_exec = config.get("samtools_exec")
    if not samtools_exec and len(sys.argv) > 2:
        samtools_exec = sys.argv[2]
    samtools_exec = samtools_exec or "samtools"
    
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
    pass1_store_conn = None
    pass1_store_path = None
    indexing_started = time.perf_counter()
    try:
        if not os.path.exists(bam_file + ".bai"):
            print(f"Indexing BAM {bam_file} for parallel processing...")
            pysam.index(bam_file)
    except Exception as e:
        print(f"BAM Indexing failed ({e}). Assuming BAM is not coordinate sorted.")
        print(f"Sorting input BAM to temporary file using {threads} threads...")
        temp_sorted_bam = bam_file + ".temp_sorted.bam"
        sort_cmd = [samtools_exec, "sort", "-@", str(threads), "-o", temp_sorted_bam, bam_file]
        try:
            subprocess.check_call(sort_cmd)
            print("Indexing temporary sorted BAM...")
            pysam.index(temp_sorted_bam)
        except Exception:
            # The outer finally only runs once the second try block is entered,
            # so clean up the partial temp BAM here before propagating.
            for leftover in (temp_sorted_bam, temp_sorted_bam + ".bai"):
                if os.path.exists(leftover):
                    os.remove(leftover)
            raise
        bam_file = temp_sorted_bam
    log_phase("prepare indexed BAM", indexing_started)

    try:
        # Get Chromosomes
        with pysam.AlignmentFile(bam_file, "rb") as b:
            references = list(b.references)
            mapped_counts = {stat.contig: stat.mapped for stat in b.get_index_statistics()}

        performance_opts = config.get("performance_opts", {}) or {}
        worker_count = resolve_worker_count(
            threads,
            len(references) or 1,
            performance_opts,
        )
        ref_chunks = balance_reference_chunks(references, mapped_counts, worker_count)
        
        print(f"Pass 1: Parallel Counting ({worker_count} workers, {len(ref_chunks)} chunks)...")
        
        # --- PASS 1: Parallel Counting ---
        # Store each chromosome-chunk result on disk. The previous in-memory
        # reduce kept every raw UMI observation alive until clustering began,
        # making the parent-process peak grow with the complete experiment.
        pass1_store_conn, pass1_store_path = open_pass1_store(
            out_dir,
            project,
            performance_opts.get("tmp_root"),
        )
        print(f"Pass 1 temporary store: {pass1_store_path}")

        pass1_args = [
            (bam_file, chunk, barcode_set, gene_set, count_introns, collect_global_umis)
            for chunk in ref_chunks
        ]
        expected_chunks = len(pass1_args)
        completed_chunks = 0
        res = partial_read = partial_umi = partial_global = None
        pass1_started = time.perf_counter()
        with multiprocessing.Pool(worker_count) as pool:
            for res in pool.imap_unordered(count_worker, pass1_args):
                completed_chunks += 1
                partial_read, partial_umi, partial_global = res
                store_pass1_result(
                    pass1_store_conn,
                    partial_read,
                    partial_umi,
                    partial_global if collect_global_umis else {},
                    commit=False,
                )
                del partial_read, partial_umi, partial_global, res

        # Do not keep the task argument list alive while the serialized store
        # is queried for clustering.
        del pass1_args
        if completed_chunks != expected_chunks:
            raise RuntimeError(
                f"DGE counting returned {completed_chunks}/{expected_chunks} chromosome chunks."
            )
        # The pass-1 store is private to this process.  Committing once avoids
        # one transaction boundary per chromosome chunk. Synchronous writes
        # are disabled for this disposable store.
        pass1_store_conn.commit()
        finalize_pass1_store(pass1_store_conn)
        log_phase(f"Pass 1 counting ({completed_chunks} chunks)", pass1_started)

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
        
        # Calculate UMI Counts (with Clustering)
        umi_workloads = pass1_barcode_workloads(
            pass1_store_conn,
            collect_global_umis,
        )
        cluster_barcodes = workload_order(umi_workloads)
        total_bcs = len(cluster_barcodes)
        cluster_workers = resolve_worker_count(threads, total_bcs or 1, performance_opts)
        
        print(
            f"Clustering UMIs for {total_bcs} barcodes using {cluster_workers} workers "
            f"(Ham_Dist={ham_dist}, correction_map={'on' if need_correction_map else 'off'}, "
            f"saturation={'on' if collect_global_umis else 'off'}, bounded scheduling)..."
        )
        
        correction_map = defaultdict(lambda: defaultdict(dict)) if need_correction_map else {}

        # Load payloads in the parent (SQLite connection owner), only when a
        # task slot is available. A slow well no longer stalls the next batch.
        max_pending = max(1, min(64, cluster_workers * 2))
        processed_bcs = 0
        clustering_started = time.perf_counter()

        def clustering_arguments():
            for bc in cluster_barcodes:
                umi_bundle, global_counts = load_pass1_umi_bundle(
                    pass1_store_conn, bc, include_global=collect_global_umis,
                )
                yield (
                    bc, umi_bundle["exon"], umi_bundle["intron"],
                    global_counts if collect_global_umis else {}, ham_dist,
                    need_correction_map, collect_global_umis, collect_global_umis,
                )

        with ProcessPoolExecutor(max_workers=cluster_workers) as pool:
            for res in bounded_results(pool, cluster_with_global, clustering_arguments(), max_pending):
                processed_bcs += 1
                if processed_bcs % 100 == 0 or processed_bcs == total_bcs:
                    print(f"Clustering {processed_bcs}/{total_bcs}...", end='\r')
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
        log_phase(f"UMI clustering ({processed_bcs} barcodes)", clustering_started)

        # Keep the store open until read matrices are rebuilt below.  Delaying
        # that pass means the large read-count maps do not coexist with the
        # UMI-clustering maps, which lowers the parent-process peak memory.
        del umi_workloads
        del cluster_barcodes

        print("\nWriting Matrices...")
        matrix_started = time.perf_counter()
        
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
            _write_json_atomic(sat_file, dict(global_umi_freq))
            _write_json_atomic(gene_sat_file, dict(gene_umi_freq))

        # Write UMI matrices first, then release their large maps before
        # reconstructing read matrices from the disk-backed pass-1 store.
        write_sparse_matrix(final_umi_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.umi")
        if count_introns:
            write_sparse_matrix(final_umi_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.umi")
            write_sparse_matrix(final_umi_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.umi")

        umi_matrix_stats = summarize_exon_intron_counts(
            final_umi_counts['exon'], final_umi_counts['intron']
        )
        del final_umi_counts

        final_read_counts = {
            'exon': defaultdict(lambda: defaultdict(int)),
            'intron': defaultdict(lambda: defaultdict(int)),
            'inex': defaultdict(lambda: defaultdict(int))
        }

        # Rebuild read matrices one barcode at a time from the same store.
        # This pass is intentionally after UMI output so both matrix families
        # are not resident in memory at the same time.
        for bc in pass1_barcodes(pass1_store_conn, "read"):
            read_bundle = load_pass1_read_bundle(pass1_store_conn, bc)
            exon_counts = read_bundle["exon"]
            intron_counts = read_bundle["intron"]
            if exon_counts:
                final_read_counts['exon'][bc].update(exon_counts)
            if intron_counts:
                final_read_counts['intron'][bc].update(intron_counts)
            for gene, count in exon_counts.items():
                final_read_counts['inex'][bc][gene] += count
            for gene, count in intron_counts.items():
                final_read_counts['inex'][bc][gene] += count

        write_sparse_matrix(final_read_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.read")
        if count_introns:
            write_sparse_matrix(final_read_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.read")
            write_sparse_matrix(final_read_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.read")

        cell_matrix_stats = {
            "schema_version": 1,
            "read_count_unit": "read_pairs" if str(config.get("read_layout", "PE")).upper() == "PE" else "reads",
            "umi": umi_matrix_stats,
            "read": summarize_exon_intron_counts(final_read_counts['exon'], final_read_counts['intron']),
        }
        cell_stats_path = os.path.join(stats_dir(out_dir), f"{project}.cell_matrix_stats.json")
        _write_json_atomic(cell_stats_path, cell_matrix_stats)
        print(f"Cell matrix QC summary written: {cell_stats_path}")
        log_phase("write expression matrices", matrix_started)

        # Matrix files and the compact QC summary are now the durable outputs.
        # Release the large in-memory count maps before optional H5AD export,
        # which reads the matrices back into sparse structures.
        del final_read_counts
        close_pass1_store(pass1_store_conn, pass1_store_path)
        pass1_store_conn = None
        pass1_store_path = None

        if _as_bool(config.get('make_h5ad', True), default=True):
            h5ad_started = time.perf_counter()
            print("Exporting combined H5AD...")
            h5ad_path = export_h5ad(out_dir, project, config=config)
            print(f"H5AD written: {h5ad_path}")
            log_phase("export H5AD", h5ad_started)

        if not make_sorted_bam and not make_ub_bam:
            log_phase("DGE total", analysis_started)
            return

        if make_ub_bam and not make_sorted_bam:
            if not out_bam:
                out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.bam")
            print(f"Writing UB-corrected BAM in a single pass ({threads} threads for BAM IO)...")
            ub_started = time.perf_counter()
            write_corrected_bam_single_pass(bam_file, out_bam, correction_map, ham_dist, threads)
            log_phase("write UB-corrected BAM", ub_started)
            log_phase("DGE total", analysis_started)
            return

        # Input BAM is coordinate-sorted here either because it already had a
        # valid index or because we created a temporary sorted copy above, so a
        # single-pass write preserves sorted order without per-chrom chunk BAMs.
        print(f"Writing UB-corrected sorted BAM in a single pass ({threads} threads for BAM IO)...")
        ub_started = time.perf_counter()
        write_corrected_bam_single_pass(bam_file, out_bam, correction_map, ham_dist, threads)

        print("Indexing Final BAM...")
        pysam.index(out_bam)
        log_phase("write/index UB-corrected BAM", ub_started)
        log_phase("DGE total", analysis_started)
        
    finally:
        if pass1_store_conn is not None or pass1_store_path:
            close_pass1_store(pass1_store_conn, pass1_store_path)
        if temp_sorted_bam and os.path.exists(temp_sorted_bam):
            print(f"Cleaning up temporary sorted BAM: {temp_sorted_bam}")
            os.remove(temp_sorted_bam)
            if os.path.exists(temp_sorted_bam + ".bai"):
                os.remove(temp_sorted_bam + ".bai")

def main():
    if pysam is None:
        raise RuntimeError("DGE analysis requires pysam. Install project dependencies first.")
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
    make_sorted_bam = _as_bool(config.get('make_sorted_bam', False))
    make_ub_bam = _as_bool(config.get('make_ub_bam', False))
    out_bam = None
    if make_sorted_bam:
        out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
    elif make_ub_bam:
        out_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.bam")
    process_bam_and_matrix(
        input_bam,
        out_bam,
        config,
        threads=num_threads,
        samtools_exec=config.get("samtools_exec"),
    )
    print("DGE Analysis pipeline finished.")

if __name__ == "__main__":
    main()
