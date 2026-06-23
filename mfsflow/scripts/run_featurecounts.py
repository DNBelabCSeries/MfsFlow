#!/usr/bin/env python3
"""
FeatureCounts execution for standalone script execution.

This module handles execution and processing of featureCounts for gene
expression quantification, including count matrix generation and
statistical analysis for single-cell RNA sequencing data.
"""

import sys
import os
import yaml
import subprocess
import csv
import shutil
import json
import collections
import gzip
import glob
import hashlib
from functools import lru_cache

try:
    from mfsflow.scripts.path_layout import expression_dir, stats_dir, load_config
except ImportError:
    from path_layout import expression_dir, stats_dir, load_config

# Constants
LOG_INTERVAL_READS = 10_000_000  # Print progress every N reads
GENEBODY_BIN_SIZE = 100_000  # Bin size for gene body coverage calculation
LRU_CACHE_SIZE = 65_536  # Default LRU cache size for parsing functions

def check_dependencies(samtools_exec, featurecounts_exec):
    def check_one(tool, name):
        if os.path.isabs(tool) or os.path.sep in str(tool):
            if not (os.path.exists(tool) and os.access(tool, os.X_OK)):
                print(f"Error: {name} is not executable: {tool}")
                sys.exit(1)
        else:
            if shutil.which(tool) is None:
                print(f"Error: {name} is not found in PATH: {tool}")
                sys.exit(1)

    check_one(samtools_exec, "samtools")
    check_one(featurecounts_exec, "featureCounts")

def get_bam_chromosomes(bam_file, samtools_exec='samtools'):
    """Reads chromosome names from BAM header."""
    cmd = [samtools_exec, 'view', '-H', bam_file]
    try:
        output = subprocess.check_output(cmd, universal_newlines=True)
        chroms = set()
        for line in output.splitlines():
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('SN:'):
                        chroms.add(part[3:])
        return chroms
    except subprocess.CalledProcessError:
        print(f"Warning: Could not read header from {bam_file}")
        return set()


def estimate_primary_mapped_reads(bam_file, samtools_exec='samtools'):
    """Count primary mapped alignments for coverage sampling fraction estimation."""
    try:
        output = subprocess.check_output(
            [samtools_exec, "view", "-c", "-F", "2308", bam_file],
            text=True,
        )
        return int(str(output).strip() or 0)
    except Exception as exc:
        print(f"Warning: Could not estimate mapped read count for {bam_file}: {exc}")
        return 0


def coverage_sampling_fraction(total_reads, target_reads):
    if not target_reads or int(target_reads) <= 0:
        return 1.0
    total_reads = int(total_reads or 0)
    if total_reads <= 0:
        return 1.0
    return min(1.0, float(target_reads) / float(total_reads))


def should_sample_for_coverage(read_name, sample_fraction, seed=42, source_label=None):
    if sample_fraction >= 1.0:
        return True
    if sample_fraction <= 0:
        return False
    key = f"{seed}|{source_label or ''}|{read_name or ''}".encode("utf-8", errors="ignore")
    value = int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big")
    return (value / float(1 << 64)) < sample_fraction

def load_gene_models(gtf_file):
    """
    Loads transcript-body models from GTF for RSeQC-like gene body coverage.
    For each gene, the longest transcript is selected, its exons are merged, and
    read overlaps are later projected onto 100 transcript-body bins in 5'->3'
    orientation.
    """
    print(f"Loading gene models for stats from {gtf_file}...")
    transcripts = collections.defaultdict(lambda: collections.defaultdict(lambda: {
        "chrom": None,
        "strand": ".",
        "exons": [],
    }))

    if not os.path.exists(gtf_file):
        print("Warning: GTF file not found. Coverage stats will be skipped.")
        return {}

    opener = gzip.open if str(gtf_file).endswith(".gz") else open
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            gene_id = _parse_gene_id(attributes)
            transcript_id = _parse_transcript_id(attributes) or gene_id

            if gene_id and transcript_id:
                tx = transcripts[gene_id][transcript_id]
                tx["chrom"] = chrom
                tx["strand"] = strand
                tx["exons"].append((start, end + 1))

    models = {}
    interval_index = collections.defaultdict(list)
    for gene_id, tx_by_id in transcripts.items():
        best = None
        for transcript_id, tx in tx_by_id.items():
            exons = _merge_half_open_intervals(tx["exons"])
            length = sum(e - s for s, e in exons)
            if length < 100:
                continue
            candidate = {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "chrom": tx["chrom"],
                "strand": tx["strand"],
                "exons": exons,
                "length": length,
            }
            if best is None or (length, transcript_id) > (best["length"], best["transcript_id"]):
                best = candidate
        if best is None:
            continue

        models[gene_id] = best
        for exon_start, exon_end in best["exons"]:
            first_bin = exon_start // _GENEBODY_BIN_SIZE
            last_bin = (exon_end - 1) // _GENEBODY_BIN_SIZE
            for bin_id in range(first_bin, last_bin + 1):
                interval_index[(best["chrom"], bin_id)].append(gene_id)

    models["__index__"] = interval_index
    print(f"Loaded {len(models) - 1} gene models.")
    return models


_GENEBODY_BIN_SIZE = GENEBODY_BIN_SIZE


@lru_cache(maxsize=LRU_CACHE_SIZE)
def _parse_transcript_id(attributes):
    result = None
    if 'transcript_id "' in attributes:
        result = attributes.split('transcript_id "')[1].split('"')[0]
    elif 'transcript_id' in attributes:
        try:
            result = attributes.split('transcript_id')[1].strip().split(';')[0].strip().strip('"')
        except Exception:
            result = None
    return result


def _merge_half_open_intervals(intervals):
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = []
    curr_start, curr_end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= curr_end:
            curr_end = max(curr_end, next_end)
        else:
            merged.append((curr_start, curr_end))
            curr_start, curr_end = next_start, next_end
    merged.append((curr_start, curr_end))
    return merged


@lru_cache(maxsize=LRU_CACHE_SIZE)
def _parse_gene_id(attributes):
    result = None
    if 'gene_id "' in attributes:
        result = attributes.split('gene_id "')[1].split('"')[0]
    elif 'gene_id' in attributes:
        try:
            result = attributes.split('gene_id')[1].strip().split(';')[0].strip().strip('"')
        except Exception:
            result = None
    return result


@lru_cache(maxsize=LRU_CACHE_SIZE)
def _parse_gene_name(attributes, gene_id):
    result = gene_id
    if 'gene_name "' in attributes:
        result = attributes.split('gene_name "')[1].split('"')[0]
    return result


def _merge_intervals(intervals):
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = []
    curr_start, curr_end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= curr_end + 1:
            curr_end = max(curr_end, next_end)
        else:
            merged.append((curr_start, curr_end))
            curr_start, curr_end = next_start, next_end
    merged.append((curr_start, curr_end))
    return merged


def _intersect_intervals(left, right):
    out = []
    i = 0
    j = 0
    left = _merge_intervals(left)
    right = _merge_intervals(right)
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if start <= end:
            out.append((start, end))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return out


def _make_global_gaps(intervals):
    merged = _merge_intervals(intervals)
    gaps = []
    for i in range(len(merged) - 1):
        start = merged[i][1] + 1
        end = merged[i + 1][0] - 1
        if start <= end:
            gaps.append((start, end))
    return gaps


def _disjoin_gene_bodies(genes):
    by_locus = collections.defaultdict(list)
    for gene_id, data in genes.items():
        if not data["exons"]:
            continue
        starts = [x[0] for x in data["exons"]]
        ends = [x[1] for x in data["exons"]]
        by_locus[(data["chrom"], data["strand"])].append((min(starts), max(ends), gene_id))

    unique_segments = collections.defaultdict(list)
    for (chrom, strand), bodies in by_locus.items():
        breakpoints = set()
        for start, end, _gene_id in bodies:
            breakpoints.add(start)
            breakpoints.add(end + 1)
        ordered = sorted(breakpoints)
        for left, right_next in zip(ordered, ordered[1:]):
            seg_start = left
            seg_end = right_next - 1
            if seg_start > seg_end:
                continue
            covering = [
                gene_id
                for body_start, body_end, gene_id in bodies
                if body_start <= seg_start and seg_end <= body_end
            ]
            if len(covering) == 1:
                unique_segments[covering[0]].append((seg_start, seg_end))
    return {gene_id: _merge_intervals(parts) for gene_id, parts in unique_segments.items()}


def _build_introns_like_r(genes):
    global_exons = collections.defaultdict(list)
    for data in genes.values():
        merged = _merge_intervals(data["exons"])
        global_exons[(data["chrom"], data["strand"])].extend(merged)

    global_gaps = {
        locus: _make_global_gaps(exons)
        for locus, exons in global_exons.items()
    }
    unique_gene_segments = _disjoin_gene_bodies(genes)

    introns = {}
    for gene_id, data in genes.items():
        locus = (data["chrom"], data["strand"])
        gaps = global_gaps.get(locus, [])
        unique_body = unique_gene_segments.get(gene_id, [])
        gene_introns = _intersect_intervals(gaps, unique_body)
        introns[gene_id] = [
            (start, end)
            for start, end in gene_introns
            if 10 < (end - start + 1) < 100000
        ]
    return introns


def parse_gtf_and_create_saf(gtf_file, out_prefix, valid_chroms=None, exon_extension=False, extension_length=0, buffer_length=0):
    """
    Parse GTF and create a combined exon/intron SAF.

    Introns are exon gaps clipped to unique gene-body regions, approximating the
    GenomicRanges/plyranges construction used by the original workflow.
    Returns path to combined_saf, and a dictionary mapping gene_id to gene_name.
    """
    print(f"Parsing GTF: {gtf_file}...")

    genes = {}
    gene_id_to_name = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue

            feature_type = parts[2]
            if feature_type != 'exon':
                continue

            chrom = parts[0]
            if valid_chroms and chrom not in valid_chroms:
                continue

            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            gene_id = _parse_gene_id(attributes)
            if not gene_id:
                continue

            gene_name = _parse_gene_name(attributes, gene_id)

            # Store mapping
            if gene_id not in gene_id_to_name:
                gene_id_to_name[gene_id] = gene_name

            if gene_id not in genes:
                genes[gene_id] = {'chrom': chrom, 'strand': strand, 'exons': []}

            genes[gene_id]['exons'].append((start, end))

    print(f"Loaded {len(genes)} genes. Generating Combined SAF file...")

    if exon_extension and extension_length > 0:
        print(f"Applying exon extension by {extension_length}bp with {buffer_length}bp buffer...")
        # Collect all merged exons globally by locus to find downstream/upstream gaps
        all_exons_by_strand = collections.defaultdict(list)
        for g_id, data in genes.items():
            for e_start, e_end in data['exons']:
                all_exons_by_strand[(data['chrom'], data['strand'])].append((e_start, e_end))

        all_merged_exons = {}
        for locus, e_list in all_exons_by_strand.items():
            all_merged_exons[locus] = _merge_intervals(e_list)

        for gene_id, data in genes.items():
            if not data['exons']: continue
            chrom = data['chrom']
            strand = data['strand']
            global_exons = all_merged_exons[(chrom, strand)]
            exons = data['exons']

            if strand == '+':
                # Find the 3' most exon (max end)
                target_idx = max(range(len(exons)), key=lambda i: exons[i][1])
                e_start, e_end = exons[target_idx]

                # Find the closest downstream exon across all genes
                next_start = None
                for (ge_start, ge_end) in global_exons:
                    if ge_start > e_end:
                        next_start = ge_start
                        break

                if next_start is not None:
                    dist = next_start - e_end
                    if dist > extension_length + buffer_length:
                        ext = extension_length
                    elif dist <= buffer_length:
                        ext = 0
                    else:
                        ext = dist - buffer_length
                    exons[target_idx] = (e_start, e_end + ext)
                else:
                    exons[target_idx] = (e_start, e_end + extension_length)

            elif strand == '-':
                # Find the 3' most exon (min start)
                target_idx = min(range(len(exons)), key=lambda i: exons[i][0])
                e_start, e_end = exons[target_idx]

                # Find the closest upstream exon across all genes (reverse sorted)
                prev_end = None
                for (ge_start, ge_end) in reversed(global_exons):
                    if ge_end < e_start:
                        prev_end = ge_end
                        break

                if prev_end is not None:
                    dist = e_start - prev_end
                    if dist > extension_length + buffer_length:
                        ext = extension_length
                    elif dist <= buffer_length:
                        ext = 0
                    else:
                        ext = dist - buffer_length
                    exons[target_idx] = (e_start - ext, e_end)
                else:
                    exons[target_idx] = (e_start - extension_length, e_end)

    combined_saf_path = f"{out_prefix}.combined.saf"
    exon_saf_path = f"{out_prefix}.exon.saf"
    intron_saf_path = f"{out_prefix}.intron.saf"
    introns_by_gene = _build_introns_like_r(genes)

    with open(combined_saf_path, 'w') as f_out, \
         open(exon_saf_path, 'w') as f_exon, \
         open(intron_saf_path, 'w') as f_intron:
        header = "GeneID\tChr\tStart\tEnd\tStrand\n"
        f_out.write(header)
        f_exon.write(header)
        f_intron.write(header)

        for gene_id, data in genes.items():
            chrom = data['chrom']
            strand = data['strand']
            merged = _merge_intervals(data['exons'])

            # Write Exons
            for start, end in merged:
                line = f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n"
                f_out.write(line)
                f_exon.write(line)

            for intron_start, intron_end in introns_by_gene.get(gene_id, []):
                line = f"{gene_id}__INTRON__\t{chrom}\t{intron_start}\t{intron_end}\t{strand}\n"
                f_out.write(line)
                f_intron.write(line)

    return combined_saf_path, gene_id_to_name


def saf_paths_from_prefix(out_prefix):
    return {
        "combined": f"{out_prefix}.combined.saf",
        "exon": f"{out_prefix}.exon.saf",
        "intron": f"{out_prefix}.intron.saf",
    }

def build_featurecounts_cmd(
    featurecounts_exec,
    input_bam,
    saf_file,
    output_counts,
    threads,
    strand_mode,
    read_layout,
    fraction_overlap=0,
    allow_multi_overlap=False,
    tmp_dir=None,
):
    layout = str(read_layout or "PE").upper()
    cmd = [
        featurecounts_exec,
        '-M',
        '-a', saf_file,
        '-F', 'SAF',
        '-o', output_counts,
        '-T', str(max(1, int(threads))),
        '-R', 'BAM',
        '-s', str(strand_mode),
        '--primary',
        '-Q', '0'
    ]
    if layout == "PE":
        cmd.extend(['-p', '-C'])
    if allow_multi_overlap:
        cmd.append('-O')
    if float(fraction_overlap) > 0:
        cmd.extend(['--fracOverlap', str(fraction_overlap)])
    if tmp_dir:
        cmd.extend(['--tmpDir', str(tmp_dir)])
    cmd.extend(['--largestOverlap', input_bam])
    return cmd

def run_featurecounts_cmd(
    featurecounts_exec,
    input_bam,
    saf_file,
    out_prefix,
    threads,
    strand_mode,
    feature_type,
    read_layout="PE",
    fraction_overlap=0,
    allow_multi_overlap=False,
    tmp_dir=None,
):
    """
    Runs featureCounts with optimized temporary file handling.
    """
    threads = max(1, int(threads))
    print(f"Running featureCounts for {feature_type} (Strand: {strand_mode}, Layout: {read_layout}, Threads: {threads})...")
    output_counts = f"{out_prefix}.counts.txt"

    cmd = build_featurecounts_cmd(
        featurecounts_exec,
        input_bam,
        saf_file,
        output_counts,
        threads,
        strand_mode,
        read_layout,
        fraction_overlap=fraction_overlap,
        allow_multi_overlap=allow_multi_overlap,
        tmp_dir=tmp_dir,
    )

    subprocess.check_call(cmd)

    for suffix in ["", ".summary"]:
        path = output_counts + suffix
        if os.path.exists(path):
            os.remove(path)

    generated_bam = f"{input_bam}.featureCounts.bam"
    if not os.path.exists(generated_bam):
        raise FileNotFoundError(f"featureCounts did not generate {generated_bam}")

    target_bam = f"{out_prefix}.bam"
    os.rename(generated_bam, target_bam)

    return target_bam


def run_featurecounts_r_order(
    featurecounts_exec,
    input_bam,
    exon_saf,
    intron_saf,
    out_prefix,
    threads,
    strand_mode,
    source_label,
    read_layout="PE",
    count_introns=True,
    fraction_overlap=0,
    allow_multi_overlap=False,
    tmp_dir=None,
):
    """
    Run featureCounts in R-order mode with optimized file handling.

    Uses a single featureCounts pass with combined SAF when possible,
    or falls back to two-pass with minimal intermediate files.
    """
    if not count_introns:
        return run_featurecounts_cmd(
            featurecounts_exec, input_bam, exon_saf, f"{out_prefix}.exon",
            threads, strand_mode, f"{source_label}_Exon",
            read_layout=read_layout, fraction_overlap=fraction_overlap,
            allow_multi_overlap=allow_multi_overlap, tmp_dir=tmp_dir,
        )

    exon_bam = run_featurecounts_cmd(
        featurecounts_exec, input_bam, exon_saf, f"{out_prefix}.exon",
        threads, strand_mode, f"{source_label}_Exon",
        read_layout=read_layout, fraction_overlap=fraction_overlap,
        allow_multi_overlap=allow_multi_overlap, tmp_dir=tmp_dir,
    )

    exon_assignments = _collect_exon_assignments(exon_bam, threads)

    intron_bam = run_featurecounts_cmd(
        featurecounts_exec, exon_bam, intron_saf, f"{out_prefix}.intron",
        threads, strand_mode, f"{source_label}_Intron",
        read_layout=read_layout, fraction_overlap=fraction_overlap,
        allow_multi_overlap=allow_multi_overlap, tmp_dir=tmp_dir,
    )

    if exon_assignments:
        _restore_exon_priority(intron_bam, exon_assignments, threads)

    _cleanup_intermediate(exon_bam)
    return intron_bam


def _collect_exon_assignments(bam_path, threads):
    """Collect read names with exon assignments from featureCounts output."""
    exon_map = {}
    try:
        import pysam
        with pysam.AlignmentFile(bam_path, "rb", threads=int(threads)) as bam:
            for read in bam:
                xt_tags = [v for t, v in read.get_tags() if t == "XT"]
                for xt_val in xt_tags:
                    if xt_val and "__INTRON__" not in xt_val:
                        exon_map[read.query_name] = xt_val
                        break
    except ImportError:
        pass
    return exon_map


def _restore_exon_priority(bam_path, exon_assignments, threads):
    """Restore exon XT tags for reads that were assigned to exons in the first pass."""
    try:
        import pysam
        tmp_bam = bam_path + ".tmp"
        with pysam.AlignmentFile(bam_path, "rb", threads=int(threads)) as infile:
            with pysam.AlignmentFile(tmp_bam, "wb", template=infile, threads=int(threads)) as outfile:
                for read in infile:
                    if read.query_name in exon_assignments:
                        read.set_tag("XT", exon_assignments[read.query_name])
                    outfile.write(read)
        os.replace(tmp_bam, bam_path)
    except ImportError:
        pass


def _cleanup_intermediate(path):
    """Remove intermediate file if it exists."""
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


def _cleanup_featurecounts_temp(input_bam):
    """Clean up featureCounts temporary files for a given input BAM."""
    if not input_bam:
        return
    patterns = [
        f"{input_bam}.featureCounts.bam",
        f"{input_bam}.featureCounts.bam.tmp",
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass


def cleanup_featurecounts_intermediates(input_bams):
    for input_bam in input_bams:
        if not input_bam:
            continue
        for path in glob.glob(f"{input_bam}.fc*"):
            if os.path.exists(path):
                os.remove(path)


_INTRON_BIN_SIZE = 100000


def load_saf_interval_index(saf_file):
    index = collections.defaultdict(list)
    if not saf_file or not os.path.exists(saf_file):
        return {}, {}
    with open(saf_file, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = str(row.get("GeneID") or "").strip()
            chrom = str(row.get("Chr") or "").strip()
            strand = str(row.get("Strand") or ".").strip() or "."
            if not gene_id or not chrom:
                continue
            if "__INTRON__" in gene_id:
                gene_id = gene_id.split("__INTRON__")[0]
            try:
                start = int(row["Start"]) - 1
                end = int(row["End"])
            except Exception:
                continue
            if start < end:
                first_bin = start // _INTRON_BIN_SIZE
                last_bin = (end - 1) // _INTRON_BIN_SIZE
                for bin_id in range(first_bin, last_bin + 1):
                    index[(chrom, strand, bin_id)].append((start, end, gene_id))

    for values in index.values():
        values.sort()
    return dict(index), {}


def resolve_counting_strand_modes(counting_opts):
    def _coerce(value, default):
        try:
            ivalue = int(value)
        except Exception:
            ivalue = default
        return ivalue if ivalue in (0, 1, 2) else default

    umi_mode = _coerce(counting_opts.get("strand", 1), 1)
    internal_mode = _coerce(counting_opts.get("internal_strand", 0), 0)
    return umi_mode, internal_mode


def should_count_read(read_obj):
    if isinstance(read_obj, str):
        parts = read_obj.split('\t')
        if len(parts) <= 1:
            return False
        flag = int(parts[1])
        return (flag & 0x900) == 0
    return not (read_obj.is_secondary or read_obj.is_supplementary)


READ_CATEGORIES = {"Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped"}


def normalize_read_category(category):
    return category if category in READ_CATEGORIES else "Other_Unassigned"


def _candidate_introns(chrom, strand_mode, is_reverse, intron_index):
    if strand_mode == 1:
        strands = ["-" if is_reverse else "+"]
    elif strand_mode == 2:
        strands = ["+" if is_reverse else "-"]
    else:
        strands = ["+", "-", "."]
    return [(chrom, strand) for strand in strands]


def _best_intron_assignment(chrom, blocks, strand_mode, is_reverse, intron_index, intron_starts):
    overlaps = collections.defaultdict(int)
    seen = set()
    for chrom_strand in _candidate_introns(chrom, strand_mode, is_reverse, intron_index):
        for block_start, block_end in blocks:
            first_bin = block_start // _INTRON_BIN_SIZE
            last_bin = (block_end - 1) // _INTRON_BIN_SIZE
            for bin_id in range(first_bin, last_bin + 1):
                key = (chrom_strand[0], chrom_strand[1], bin_id)
                for intron_start, intron_end, gene_id in intron_index.get(key, ()):
                    interval_key = (chrom_strand[0], chrom_strand[1], intron_start, intron_end, gene_id, block_start, block_end)
                    if interval_key in seen:
                        continue
                    seen.add(interval_key)
                    ov = min(block_end, intron_end) - max(block_start, intron_start)
                    if ov > 0:
                        overlaps[gene_id] += ov
    if not overlaps:
        return None, "Intergenic"
    ranked = sorted(overlaps.items(), key=lambda item: (-item[1], item[0]))
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None, "Ambiguity"
    return ranked[0][0], "Intron"


def _sam_blocks(pos_1based, cigar):
    blocks = []
    if not cigar or cigar == "*":
        return blocks
    curr_pos = int(pos_1based) - 1
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
            continue
        if ch in "M=X":
            blocks.append((curr_pos, curr_pos + num))
            curr_pos += num
        elif ch in "DN":
            curr_pos += num
        elif ch in "ISH":
            pass
        num = 0
    return blocks


def _pysam_blocks_1based_half_open(read_obj):
    """Convert pysam 0-based half-open blocks to 1-based half-open blocks."""
    return [(start + 1, end + 1) for start, end in read_obj.get_blocks()]


def _sam_line_blocks_1based_half_open(parts):
    pos = int(parts[3])
    cigar = parts[5]
    blocks = []
    curr_pos = pos
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
        else:
            if ch in "M=X":
                blocks.append((curr_pos, curr_pos + num))
                curr_pos += num
            elif ch in "DN":
                curr_pos += num
            elif ch in "SH":
                pass
            num = 0
    return blocks


def _candidate_gene_body_models(chrom, blocks, gene_models):
    index = gene_models.get("__index__", {})
    candidate_ids = set()
    for block_start, block_end in blocks:
        first_bin = block_start // _GENEBODY_BIN_SIZE
        last_bin = (block_end - 1) // _GENEBODY_BIN_SIZE
        for bin_id in range(first_bin, last_bin + 1):
            candidate_ids.update(index.get((chrom, bin_id), ()))
    return [gene_models[gene_id] for gene_id in candidate_ids if gene_id in gene_models]


def _project_blocks_to_gene_body(model, blocks):
    increments = [0.0] * 100
    total_overlap = 0
    length = model["length"]

    exon_offsets = []
    offset = 0
    exons_5_to_3 = list(reversed(model["exons"])) if model["strand"] == "-" else model["exons"]
    for exon_start, exon_end in exons_5_to_3:
        exon_offsets.append((exon_start, exon_end, offset))
        offset += exon_end - exon_start

    for block_start, block_end in blocks:
        for exon_start, exon_end, exon_offset in exon_offsets:
            ov_start = max(block_start, exon_start)
            ov_end = min(block_end, exon_end)
            if ov_end <= ov_start:
                continue

            if model["strand"] == "-":
                tx_start = exon_offset + (exon_end - ov_end)
                tx_end = exon_offset + (exon_end - ov_start)
            else:
                tx_start = exon_offset + (ov_start - exon_start)
                tx_end = exon_offset + (ov_end - exon_start)

            total_overlap += ov_end - ov_start
            _add_transcript_interval_to_bins(tx_start, tx_end, length, increments)

    return total_overlap, increments


def _add_transcript_interval_to_bins(tx_start, tx_end, length, increments):
    if tx_end <= tx_start:
        return
    first_bin = min(99, (tx_start * 100) // length)
    last_bin = min(99, ((tx_end - 1) * 100) // length)
    for bin_idx in range(first_bin, last_bin + 1):
        bin_start = (bin_idx * length) // 100
        bin_end = ((bin_idx + 1) * length) // 100
        if bin_end <= bin_start:
            continue
        overlap = min(tx_end, bin_end) - max(tx_start, bin_start)
        if overlap > 0:
            increments[bin_idx] += overlap / float(bin_end - bin_start)

def process_bam_and_calculate_stats(
    input_bam,
    out_bam,
    samtools_exec,
    threads=4,
    gene_map=None,
    source_label=None,
    gene_models=None,
    collect_coverage=True,
    output_handle=None,
    intron_index=None,
    intron_starts=None,
    strand_mode=0,
    coverage_sample_fraction=1.0,
    coverage_sample_seed=42,
):
    """
    Processes a single BAM from Combined SAF featureCounts.
    Parses XT tag to distinguish Exon (GeneID) vs Intron (GeneID__INTRON__).
    Adds RE:Z:E/N/I tag.
    Adds GN:Z:GeneName tag.
    Calculates Stats on the fly.
    """
    target_desc = "shared output handle" if output_handle is not None else out_bam
    print(f"Processing BAM {input_bam} -> {target_desc} (Source: {source_label}, Threads: {threads})...")

    read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_arr = [0] * 100
    cov_count = 0
    intron_index = intron_index or {}
    intron_starts = intron_starts or {}

    def update_stats(read_obj, category, source_lbl):
        category = normalize_read_category(category)
        bc = None
        if isinstance(read_obj, str):
            if not should_count_read(read_obj):
                return
            cb_idx = read_obj.find("CB:Z:")
            if cb_idx != -1:
                end_cb = read_obj.find('\t', cb_idx)
                if end_cb == -1: end_cb = len(read_obj)
                bc = read_obj[cb_idx+5:end_cb].strip()

            if bc:
                read_stats[bc][category] += 1
                if source_lbl == 'UMI': read_stats[bc]['UMI_Reads'] += 1
                elif source_lbl == 'Internal': read_stats[bc]['Internal_Reads'] += 1

        else: # pysam
            if not should_count_read(read_obj):
                return
            if read_obj.has_tag("CB"):
                bc = read_obj.get_tag("CB")
                read_stats[bc][category] += 1
                if source_lbl == 'UMI': read_stats[bc]['UMI_Reads'] += 1
                elif source_lbl == 'Internal': read_stats[bc]['Internal_Reads'] += 1
            else:
                 read_stats["__NO_CB__"]["Unused BC"] += 1

    def update_coverage(read_obj, gene_models, cov_arr):
        if not collect_coverage or not gene_models: return False
        if isinstance(read_obj, str):
            if not should_count_read(read_obj):
                return False
            parts = read_obj.split('\t')
            if len(parts) < 6 or parts[2] == "*" or parts[5] == "*":
                return False
            if not should_sample_for_coverage(parts[0], coverage_sample_fraction, coverage_sample_seed, source_label):
                return False
            chrom = parts[2]
            blocks = _sam_line_blocks_1based_half_open(parts)
        else:
            if not should_count_read(read_obj) or read_obj.is_unmapped:
                return False
            if not should_sample_for_coverage(read_obj.query_name, coverage_sample_fraction, coverage_sample_seed, source_label):
                return False
            chrom = read_obj.reference_name
            blocks = _pysam_blocks_1based_half_open(read_obj)

        if not blocks: return False

        candidates = _candidate_gene_body_models(chrom, blocks, gene_models)
        if not candidates:
            return False

        ranked = []
        for model in candidates:
            overlap, increments = _project_blocks_to_gene_body(model, blocks)
            if overlap > 0:
                ranked.append((overlap, model["gene_id"], increments))
        if not ranked:
            return False

        ranked.sort(key=lambda item: (-item[0], item[1]))
        if len(ranked) > 1 and ranked[0][0] == ranked[1][0]:
            return False

        for i, value in enumerate(ranked[0][2]):
            cov_arr[i] += value
        return True

    def choose_assignment_from_xt(xt_values):
        exon_gene = None
        intron_gene = None
        for xt_val in xt_values:
            if not xt_val:
                continue
            if "__INTRON__" in xt_val:
                if intron_gene is None:
                    intron_gene = xt_val.split("__INTRON__")[0]
            elif exon_gene is None:
                exon_gene = xt_val
        if exon_gene:
            return exon_gene, "E", "Exon"
        if intron_gene:
            return intron_gene, "N", "Intron"
        return None, None, None

    # Try pysam
    try:
        import pysam
        print("Using pysam for BAM processing...")

        # Input is likely name sorted from featureCounts
        with pysam.AlignmentFile(input_bam, "rb", threads=int(threads)) as f_in:
            if output_handle is None:
                f_out_ctx = pysam.AlignmentFile(out_bam, "wb", template=f_in, threads=int(threads))
            else:
                f_out_ctx = None
            f_out = output_handle if output_handle is not None else f_out_ctx
            try:
                count = 0
                for read in f_in:
                    count += 1
                    if count % LOG_INTERVAL_READS == 0: print(f"Processed {count} reads...", end='\r')

                    category = "Intergenic"
                    final_read = read

                    xt_values = [value for tag, value in read.get_tags() if tag == "XT"]
                    gene_id, re_tag, assigned_category = choose_assignment_from_xt(xt_values)
                    if gene_id:
                        final_read.set_tag('GX', gene_id)
                        final_read.set_tag('RE', re_tag)
                        category = assigned_category

                    else:
                        # Unassigned
                        status = "Unassigned"
                        if read.has_tag('XS'):
                            xs_val = read.get_tag('XS')
                            if "Unassigned_" in xs_val:
                                status = xs_val.replace("Unassigned_", "")
                            elif xs_val == "Assigned":
                                 pass

                        if status in {"NoFeatures", "Ambiguity"}:
                             intron_gene, intron_category = _best_intron_assignment(
                                 final_read.reference_name,
                                 final_read.get_blocks(),
                                 strand_mode,
                                 final_read.is_reverse,
                                 intron_index,
                                 intron_starts,
                             ) if intron_index and not final_read.is_unmapped else (None, "Intergenic")
                             if intron_gene:
                                 final_read.set_tag('GX', intron_gene)
                                 final_read.set_tag('RE', 'N')
                                 category = "Intron"
                             elif intron_category == "Ambiguity":
                                 category = "Ambiguity"
                             else:
                                 if status == "Ambiguity":
                                     category = "Ambiguity"
                                 else:
                                     final_read.set_tag('RE', 'I')
                                     category = "Intergenic"
                        elif read.is_unmapped:
                             category = "Unmapped"
                        else:
                             category = normalize_read_category(status)

                    if source_label:
                        final_read.set_tag('SR', source_label)

                    if final_read.has_tag('GX'):
                        g_id = final_read.get_tag('GX')
                        if gene_map:
                            g_name = gene_map.get(g_id, g_id)
                            final_read.set_tag('GN', g_name)

                    # Stats
                    update_stats(final_read, category, source_label)
                    if collect_coverage and update_coverage(final_read, gene_models, cov_arr):
                        cov_count += 1

                    f_out.write(final_read)
            finally:
                if f_out_ctx is not None:
                    f_out_ctx.close()

        print("\nProcessing complete (via pysam).")
        return read_stats, cov_arr, cov_count

    except ImportError:
        print("pysam not found. Falling back to samtools pipe method...")
    except Exception as e:
        print(f"pysam processing failed: {e}. Falling back to samtools pipe...")

    # Fallback to samtools pipe
    buf_size = 64 * 1024 * 1024

    cmd_in = [samtools_exec, 'view', '-h', '-@', str(threads), input_bam]
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)

    if output_handle is not None:
        raise RuntimeError("Shared BAM output handle requires pysam support.")

    cmd_out = [samtools_exec, 'view', '-b', '-@', str(threads), '-o', out_bam, '-']
    proc_out = subprocess.Popen(cmd_out, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=buf_size)

    try:
        count = 0
        for line in proc_in.stdout:
            if line.startswith('@'):
                proc_out.stdin.write(line)
                continue

            count += 1
            if count % LOG_INTERVAL_READS == 0: print(f"Processed {count} reads...", end='\r')

            final_line = None
            category = "Intergenic"

            fields = line.rstrip().split('\t')
            xt_values = [field[5:] for field in fields[11:] if field.startswith('XT:Z:')]
            gene_id, re_tag, assigned_category = choose_assignment_from_xt(xt_values)

            if gene_id:
                final_line = line.rstrip() + f"\tGX:Z:{gene_id}\tRE:Z:{re_tag}"
                category = assigned_category
            else:
                # Unassigned
                reason = None
                xs_pos = line.find('XS:Z:')
                if xs_pos != -1:
                    xs_end = line.find('\t', xs_pos)
                    reason = line[xs_pos+5 : xs_end] if xs_end != -1 else line[xs_pos+5:].rstrip()

                if reason and "Unassigned_" in reason:
                    status = reason.replace("Unassigned_", "")
                    if status in {"NoFeatures", "Ambiguity"}:
                        parts = line.rstrip().split('\t')
                        flag = int(parts[1])
                        blocks = _sam_blocks(int(parts[3]), parts[5])
                        intron_gene, intron_category = _best_intron_assignment(
                            parts[2],
                            blocks,
                            strand_mode,
                            (flag & 0x10) != 0,
                            intron_index,
                            intron_starts,
                        ) if intron_index and not (flag & 0x4) else (None, "Intergenic")
                        if intron_gene:
                            final_line = line.rstrip() + f"\tGX:Z:{intron_gene}\tRE:Z:N"
                            category = "Intron"
                        elif intron_category == "Ambiguity":
                            final_line = line.rstrip()
                            category = "Ambiguity"
                        else:
                            if status == "Ambiguity":
                                final_line = line.rstrip()
                                category = "Ambiguity"
                            else:
                                final_line = line.rstrip() + "\tRE:Z:I"
                                category = "Intergenic"
                    else:
                        final_line = line.rstrip()
                        category = normalize_read_category(status)
                else:
                    # Check unmapped flag
                    parts = line.split('\t')
                    flag = int(parts[1])
                    if not (flag & 0x4):
                         final_line = line.rstrip() + "\tRE:Z:I"
                         category = "Intergenic"
                    else:
                         final_line = line.rstrip()
                         category = "Unmapped"

            if source_label:
                final_line += f"\tSR:Z:{source_label}"

            # GN Tag
            if gene_map and ('GX:Z:' in final_line):
                start_gx = final_line.find('GX:Z:') + 5
                end_gx = final_line.find('\t', start_gx)
                gene_id = final_line[start_gx:] if end_gx == -1 else final_line[start_gx:end_gx]

                gene_name = gene_map.get(gene_id, gene_id)
                final_line += f"\tGN:Z:{gene_name}"

            # Stats
            update_stats(final_line, category, source_label)
            if collect_coverage and update_coverage(final_line, gene_models, cov_arr):
                 cov_count += 1

            proc_out.stdin.write(final_line + "\n")

    except BrokenPipeError:
        outs, errs = proc_out.communicate()
        if errs: print(errs)
        raise
    finally:
        if proc_in.stdout: proc_in.stdout.close()
        if proc_out.stdin: proc_out.stdin.close()
        proc_in.wait()
        proc_out.wait()

    print("\nProcessing complete.")
    return read_stats, cov_arr, cov_count

def split_bam_smartseq3(bam_file, threads, samtools_exec):
    print("Splitting BAM for Smart-seq3 processing (One-pass Optimized)...")
    prefix = bam_file.replace('.bam', '')
    umi_bam = f"{prefix}.UMI.bam"
    internal_bam = f"{prefix}.internal.bam"

    # Method 1: Try pysam (Fastest/Cleanest if installed)
    try:
        import pysam
        print("Using pysam for splitting...")
        read_threads = max(1, int(threads) // 2)
        write_threads = max(1, int(threads) // 4)

        with pysam.AlignmentFile(bam_file, "rb", threads=read_threads) as infile:
            with pysam.AlignmentFile(umi_bam, "wb", template=infile, threads=write_threads) as out_umi, \
                 pysam.AlignmentFile(internal_bam, "wb", template=infile, threads=write_threads) as out_int:

                for read in infile:
                    if read.has_tag('UR'):
                        val = read.get_tag('UR')
                        if val:
                            out_umi.write(read)
                        else:
                            out_int.write(read)
                    else:
                        out_int.write(read)
        return internal_bam, umi_bam
    except ImportError:
        pass # Fallback to subprocess

    # Method 2: Subprocess Pipe
    print("pysam not found, using samtools pipe...")
    buf_size = 64 * 1024 * 1024

    cmd_in = [samtools_exec, 'view', '-h', '-@', str(max(1, int(threads)//2)), bam_file]
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)

    cmd_umi = [samtools_exec, 'view', '-b', '-@', str(max(1, int(threads)//4)), '-o', umi_bam, '-']
    proc_umi = subprocess.Popen(cmd_umi, stdin=subprocess.PIPE, text=True, bufsize=buf_size)

    cmd_int = [samtools_exec, 'view', '-b', '-@', str(max(1, int(threads)//4)), '-o', internal_bam, '-']
    proc_int = subprocess.Popen(cmd_int, stdin=subprocess.PIPE, text=True, bufsize=buf_size)

    try:
        for line in proc_in.stdout:
            if line.startswith('@'):
                proc_umi.stdin.write(line)
                proc_int.stdin.write(line)
                continue

            if 'UR:Z:' in line and 'UR:Z:\t' not in line:
                proc_umi.stdin.write(line)
            else:
                proc_int.stdin.write(line)

    except BrokenPipeError:
        print("Error: Broken Pipe during BAM splitting.")
        raise
    finally:
        if proc_in.stdout: proc_in.stdout.close()
        if proc_umi.stdin: proc_umi.stdin.close()
        if proc_int.stdin: proc_int.stdin.close()
        proc_in.wait()
        proc_umi.wait()
        proc_int.wait()

    return internal_bam, umi_bam

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file')
    parser.add_argument('--umi_bam', required=False, help="Aligned UMI BAM")
    parser.add_argument('--internal_bam', required=False, help="Aligned Internal BAM")
    args = parser.parse_args()

    yaml_file = args.yaml_file
    config = load_config(yaml_file)

    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 4))
    samtools_exec = config.get('samtools_exec', 'samtools')
    featurecounts_exec = config.get('featureCounts_exec', 'featureCounts')

    check_dependencies(samtools_exec, featurecounts_exec)

    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    final_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")

    counting_opts = config.get('counting_opts', {})

    print(f"Processing Project: {project} with {num_threads} threads.")

    umi_bam = args.umi_bam
    internal_bam = args.internal_bam

    if not umi_bam or not internal_bam:
        umi_bam = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
        internal_bam = os.path.join(out_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")

    # Check existence
    if not os.path.exists(umi_bam) and not os.path.exists(internal_bam):
         print(f"Error: Input BAMs not found.")
         print(f"Checked UMI BAM: {umi_bam}")
         print(f"Checked Internal BAM: {internal_bam}")
         sys.exit(1)

    header_bam = umi_bam if os.path.exists(umi_bam) else internal_bam
    valid_chroms = get_bam_chromosomes(header_bam, samtools_exec)

    saf_dir = expression_dir(out_dir)
    if not os.path.exists(saf_dir): os.makedirs(saf_dir)

    saf_prefix = os.path.join(saf_dir, f"{project}")

    reference_opts = config.get("reference", {})
    exon_extension_flag = bool(reference_opts.get("exon_extension", False))
    extension_length = int(reference_opts.get("extension_length", 0))
    buffer_length = int(reference_opts.get("buffer_length", extension_length // 2))

    combined_saf, gene_map = parse_gtf_and_create_saf(
        gtf_file, saf_prefix, valid_chroms,
        exon_extension=exon_extension_flag,
        extension_length=extension_length,
        buffer_length=buffer_length
    )
    saf_paths = saf_paths_from_prefix(saf_prefix)

    collect_coverage = bool(config.get('make_stats', True))
    gene_models = load_gene_models(gtf_file) if collect_coverage else {}

    fc_outputs = []
    total_read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    total_cov_umi = [0] * 100
    total_cov_int = [0] * 100
    total_cov_umi_reads = 0
    total_cov_int_reads = 0
    total_cov_umi_fraction = 1.0
    total_cov_int_fraction = 1.0

    def merge_stats(dest_stats, src_stats):
        for bc, counts in src_stats.items():
            for cat, val in counts.items():
                dest_stats[bc][cat] += val

    def merge_coverage(dest_cov, src_cov):
        for i in range(100):
            dest_cov[i] += src_cov[i]

    count_introns = bool(counting_opts.get("introns", True))
    featurecounts_strategy = str(counting_opts.get("featurecounts_strategy", "hybrid")).strip().lower()
    read_layout = str(config.get('read_layout', 'PE') or 'PE').upper()
    use_r_order = featurecounts_strategy in {"r", "r_order", "two_pass", "exact"}
    fraction_overlap = counting_opts.get("fraction_overlap", 0)
    allow_multi_overlap = bool(counting_opts.get("multi_overlap", False))
    gene_body_max_reads = int(counting_opts.get("gene_body_max_reads", 50000000) or 0)
    gene_body_sample_seed = int(counting_opts.get("gene_body_sample_seed", 42) or 42)
    umi_strand_mode, internal_strand_mode = resolve_counting_strand_modes(counting_opts)
    print(
        f"featureCounts opts: strategy={featurecounts_strategy}, "
        f"multi_overlap={allow_multi_overlap}, fraction_overlap={fraction_overlap}, "
        f"read_layout={read_layout}, umi_strand={umi_strand_mode}, internal_strand={internal_strand_mode}, "
        f"gene_body_max_reads={gene_body_max_reads}, gene_body_sample_seed={gene_body_sample_seed}"
    )
    intron_index, intron_starts = ({}, {})
    if count_introns and not use_r_order:
        intron_index, intron_starts = load_saf_interval_index(saf_paths["intron"])
        print(
            f"Using hybrid counting: exon featureCounts + Python intron assignment "
            f"({sum(len(v) for v in intron_index.values())} intron intervals)."
        )
    elif use_r_order:
        print("Using R-order counting: exon featureCounts followed by intron featureCounts.")

    # Process Internal
    if os.path.exists(internal_bam):
        print("Running featureCounts for Internal Reads...")
        fc_prefix_int = internal_bam + ".fc"
        if use_r_order:
            fc_out_int = run_featurecounts_r_order(
                featurecounts_exec,
                internal_bam,
                saf_paths["exon"],
                saf_paths["intron"],
                fc_prefix_int,
                num_threads,
                internal_strand_mode,
                "Internal",
                read_layout=read_layout,
                count_introns=count_introns,
                fraction_overlap=fraction_overlap,
                allow_multi_overlap=allow_multi_overlap,
            )
        else:
            fc_out_int = run_featurecounts_cmd(
                featurecounts_exec,
                internal_bam,
                saf_paths["exon"],
                f"{fc_prefix_int}.exon",
                num_threads,
                internal_strand_mode,
                "Internal_Exon",
                read_layout=read_layout,
                fraction_overlap=fraction_overlap,
                allow_multi_overlap=allow_multi_overlap,
            )
        fc_outputs.append(("Internal", fc_out_int, internal_strand_mode))

    # Process UMI
    if os.path.exists(umi_bam):
        print("Running featureCounts for UMI Reads...")
        fc_prefix_umi = umi_bam + ".fc"
        if use_r_order:
            fc_out_umi = run_featurecounts_r_order(
                featurecounts_exec,
                umi_bam,
                saf_paths["exon"],
                saf_paths["intron"],
                fc_prefix_umi,
                num_threads,
                umi_strand_mode,
                "UMI",
                read_layout=read_layout,
                count_introns=count_introns,
                fraction_overlap=fraction_overlap,
                allow_multi_overlap=allow_multi_overlap,
            )
        else:
            fc_out_umi = run_featurecounts_cmd(
                featurecounts_exec,
                umi_bam,
                saf_paths["exon"],
                f"{fc_prefix_umi}.exon",
                num_threads,
                umi_strand_mode,
                "UMI_Exon",
                read_layout=read_layout,
                fraction_overlap=fraction_overlap,
                allow_multi_overlap=allow_multi_overlap,
            )
        fc_outputs.append(("UMI", fc_out_umi, umi_strand_mode))

    if not fc_outputs:
        print("Error: No BAMs processed.")
        sys.exit(1)

    try:
        import pysam
    except ImportError:
        pysam = None

    if pysam is not None:
        print(f"Writing final GeneTagged BAM directly: {final_bam}")
        with pysam.AlignmentFile(fc_outputs[0][1], "rb", threads=int(num_threads)) as template_in:
            with pysam.AlignmentFile(final_bam, "wb", template=template_in, threads=int(num_threads)) as final_out:
                for source_label, fc_bam, fc_strand_mode in fc_outputs:
                    mapped_for_cov = estimate_primary_mapped_reads(fc_bam, samtools_exec) if collect_coverage else 0
                    cov_fraction = coverage_sampling_fraction(mapped_for_cov, gene_body_max_reads)
                    print(
                        f"Gene body coverage sampling ({source_label}): "
                        f"mapped={mapped_for_cov}, fraction={cov_fraction:.6f}"
                    )
                    r_stats, cov, cov_reads = process_bam_and_calculate_stats(
                        fc_bam, final_bam, samtools_exec, num_threads,
                        gene_map, source_label=source_label, gene_models=gene_models,
                        collect_coverage=collect_coverage, output_handle=final_out,
                        intron_index=intron_index, intron_starts=intron_starts,
                        strand_mode=fc_strand_mode,
                        coverage_sample_fraction=cov_fraction,
                        coverage_sample_seed=gene_body_sample_seed,
                    )
                    merge_stats(total_read_stats, r_stats)
                    if source_label == "UMI":
                        merge_coverage(total_cov_umi, cov)
                        total_cov_umi_reads += cov_reads
                        total_cov_umi_fraction = cov_fraction
                    else:
                        merge_coverage(total_cov_int, cov)
                        total_cov_int_reads += cov_reads
                        total_cov_int_fraction = cov_fraction
                    os.remove(fc_bam)
    else:
        print("pysam unavailable for direct final BAM writing; falling back to intermediate BAM merge.")
        bams_to_merge = []
        for source_label, fc_bam, fc_strand_mode in fc_outputs:
            processed_bam = fc_bam + ".processed.bam"
            mapped_for_cov = estimate_primary_mapped_reads(fc_bam, samtools_exec) if collect_coverage else 0
            cov_fraction = coverage_sampling_fraction(mapped_for_cov, gene_body_max_reads)
            print(
                f"Gene body coverage sampling ({source_label}): "
                f"mapped={mapped_for_cov}, fraction={cov_fraction:.6f}"
            )
            r_stats, cov, cov_reads = process_bam_and_calculate_stats(
                fc_bam, processed_bam, samtools_exec, num_threads,
                gene_map, source_label=source_label, gene_models=gene_models,
                collect_coverage=collect_coverage,
                intron_index=intron_index, intron_starts=intron_starts,
                strand_mode=fc_strand_mode,
                coverage_sample_fraction=cov_fraction,
                coverage_sample_seed=gene_body_sample_seed,
            )
            merge_stats(total_read_stats, r_stats)
            if source_label == "UMI":
                merge_coverage(total_cov_umi, cov)
                total_cov_umi_reads += cov_reads
                total_cov_umi_fraction = cov_fraction
            else:
                merge_coverage(total_cov_int, cov)
                total_cov_int_reads += cov_reads
                total_cov_int_fraction = cov_fraction
            os.remove(fc_bam)
            bams_to_merge.append(processed_bam)

        if len(bams_to_merge) == 1:
            os.rename(bams_to_merge[0], final_bam)
        else:
            print(f"Merging BAMs with {num_threads} threads...")
            cmd = [samtools_exec, 'cat', '-@', str(num_threads), '-o', final_bam] + bams_to_merge
            subprocess.check_call(cmd)
            for b in bams_to_merge:
                if os.path.exists(b):
                    os.remove(b)

    cleanup_featurecounts_intermediates([umi_bam, internal_bam])

    # Save Stats
    stats_out = os.path.join(stats_dir(out_dir), f"{project}.read_stats.json")
    if not os.path.exists(os.path.dirname(stats_out)):
        os.makedirs(os.path.dirname(stats_out))

    stats_data = {
        "read_stats": total_read_stats,
        "coverage_umi": total_cov_umi,
        "coverage_int": total_cov_int,
        "coverage_umi_reads": total_cov_umi_reads,
        "coverage_internal_reads": total_cov_int_reads,
        "coverage_umi_sample_fraction": total_cov_umi_fraction,
        "coverage_internal_sample_fraction": total_cov_int_fraction,
        "coverage_sample_seed": gene_body_sample_seed,
        "coverage_target_reads": gene_body_max_reads
    }
    with open(stats_out, 'w') as f:
        json.dump(stats_data, f)

    print("FeatureCounts pipeline finished successfully.")

if __name__ == "__main__":
    main()
