#!/usr/bin/env python3
"""
Statistics generation for standalone script execution.

This module generates various statistics and quality metrics from processed
single-cell RNA sequencing data, including cell counts, UMI distributions,
and gene expression statistics.
"""

import sys
import os
import collections
import gzip
import json
import numpy as np

try:
    from mfsflow.path_layout import config_dir, expression_dir, stats_dir, load_config
except ImportError:
    from path_layout import config_dir, expression_dir, stats_dir, load_config

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def calculate_read_ratios(exon_reads, intron_reads, intergenic_reads, ambiguity_reads,
                          unmapped_reads, other_unassigned_reads, all_reads):
    """Return mapping, genic, and legacy exon/intron ratios."""
    # featureCounts unassigned reasons other than Unmapped still describe
    # aligned fragments (for example Chimera or FragmentLength).
    mapped_reads = (
        exon_reads
        + intron_reads
        + intergenic_reads
        + ambiguity_reads
        + other_unassigned_reads
    )
    mapping_denom = mapped_reads + unmapped_reads
    mapping_ratio = (mapped_reads / mapping_denom) if mapping_denom > 0 else None
    genic_ratio = ((exon_reads + intron_reads) / all_reads) if all_reads > 0 else None
    exon_intron_ratio = (exon_reads / intron_reads) if intron_reads > 0 else None
    return mapping_ratio, genic_ratio, exon_intron_ratio

def load_barcode_mapping(out_dir, _project):
    mapping = {}
    expect_candidates = [
        os.path.join(config_dir(out_dir), "expect_id_barcode.tsv"),
        os.path.join(out_dir, "expect_id_barcode.tsv"),
    ]
    expect_file = next((p for p in expect_candidates if os.path.exists(p)), None)
    
    if not expect_file:
        return mapping

    with open(expect_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                if parts[0] in ['wellID', 'WellID']: continue
                well = parts[0].strip()
                umi_bc = parts[1].strip()
                int_bc = parts[2].strip()
                mapping[well] = {'umi': umi_bc, 'int': int_bc}
    return mapping

def calculate_matrix_stats(matrix_dir):
    """
    Calculates UMI count and Gene count per cell from a sparse matrix directory.
    Returns: dict {bc: {'umis': int, 'genes': int, 'gene_set': set}}
    """
    stats = collections.defaultdict(lambda: {'umis': 0, 'genes': 0, 'gene_set': set()})
    
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    
    if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path)):
         matrix_path = matrix_path[:-3]
         barcodes_path = barcodes_path[:-3]
         if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path)):
             return stats

    # Load Barcodes
    opener = gzip.open if barcodes_path.endswith('.gz') else open
    with opener(barcodes_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
        
    # Parse MTX
    opener = gzip.open if matrix_path.endswith('.gz') else open
    with opener(matrix_path, 'rt') as f:
        header_passed = False
        for line in f:
            if line.startswith('%'): continue
            if not header_passed:
                header_passed = True
                continue
            
            parts = line.split()
            if len(parts) < 3: continue
            
            col_idx = int(parts[1]) - 1
            row_idx = int(parts[0]) - 1 
            val = float(parts[2])
            
            if val > 0:
                bc = barcodes[col_idx]
                stats[bc]['umis'] += int(val)
                stats[bc]['genes'] += 1
                stats[bc]['gene_set'].add(row_idx)
                
    return stats


def load_cell_matrix_stats(path):
    """Load the compact per-cell summary emitted by DGE, if available."""
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if payload.get("schema_version") != 1:
            return None
        if not isinstance(payload.get("umi"), dict) or not isinstance(payload.get("read"), dict):
            return None
        return payload
    except (OSError, ValueError, TypeError):
        return None

def plot_coverage(cov_umi, cov_int, out_prefix):
    if not HAS_MATPLOTLIB: return

    def normalize_coverage(values):
        values = np.asarray(values, dtype=np.float64).reshape(-1)
        if values.size == 100:
            return values
        normalized = np.zeros(100, dtype=np.float64)
        normalized[:min(values.size, 100)] = values[:100]
        return normalized

    cov_umi = normalize_coverage(cov_umi)
    cov_int = normalize_coverage(cov_int)
    
    def smooth(y, window=7):
        y = np.asarray(y, dtype=np.float64)
        if window <= 1: return y
        if window % 2 == 0: window += 1
        if y.size < window: return y
        kernel = np.ones(window, dtype=np.float64) / window
        pad = window // 2
        y_pad = np.pad(y, (pad, pad), mode='edge')
        return np.convolve(y_pad, kernel, mode='valid')

    sum_umi = float(np.sum(cov_umi))
    sum_int = float(np.sum(cov_int))
    cov_all = cov_umi + cov_int
    sum_all = float(np.sum(cov_all))
    # Avoid div by zero
    frac_umi = (cov_umi / sum_umi) if sum_umi > 0 else np.zeros_like(cov_umi, dtype=np.float64)
    frac_int = (cov_int / sum_int) if sum_int > 0 else np.zeros_like(cov_int, dtype=np.float64)
    frac_all = (cov_all / sum_all) if sum_all > 0 else np.zeros_like(cov_all, dtype=np.float64)

    max_umi = float(np.max(cov_umi)) if cov_umi.size else 0.0
    max_int = float(np.max(cov_int)) if cov_int.size else 0.0
    max_all = float(np.max(cov_all)) if cov_all.size else 0.0
    maxnorm_umi = (cov_umi / max_umi) if max_umi > 0 else np.zeros_like(cov_umi, dtype=np.float64)
    maxnorm_int = (cov_int / max_int) if max_int > 0 else np.zeros_like(cov_int, dtype=np.float64)
    maxnorm_all = (cov_all / max_all) if max_all > 0 else np.zeros_like(cov_all, dtype=np.float64)

    maxnorm_umi_s = smooth(maxnorm_umi, window=7)
    maxnorm_int_s = smooth(maxnorm_int, window=7)

    x = np.arange(1, 101)

    plt.figure(figsize=(10, 8))
    plt.plot(x, maxnorm_int_s, label='Internal Coverage', color='blue', linewidth=2)
    plt.plot(x, maxnorm_umi_s, label='UMI Coverage', color='red', linewidth=2)

    plt.ylim(0, 1.05)
    plt.xlabel("Gene Body Percentile (5'->3')")
    plt.ylabel("Coverage")
    plt.title("Gene Body Coverage Comparison")
    plt.grid(True, linestyle='-', alpha=0.2)
    plt.legend(frameon=False, loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.geneBodyCoverage.pdf")
    plt.close()

    with open(f"{out_prefix}.geneBodyCoverage.txt", 'w') as f:
        f.write("Percentile\tUMI_Coverage\tInternal_Coverage\tAll_Coverage\tUMI_Frac\tInternal_Frac\tAll_Frac\tUMI_MaxNorm\tInternal_MaxNorm\tAll_MaxNorm\n")
        for i in range(100):
            f.write(
                f"{i+1}\t{cov_umi[i]:.6f}\t{cov_int[i]:.6f}\t{cov_all[i]:.6f}\t"
                f"{frac_umi[i]:.6f}\t{frac_int[i]:.6f}\t{frac_all[i]:.6f}\t"
                f"{maxnorm_umi[i]:.6f}\t{maxnorm_int[i]:.6f}\t{maxnorm_all[i]:.6f}\n"
            )

def plot_gene_umi_counts_by_type(stats_exon, stats_intron, stats_inex, wells, out_pdf):
    if not HAS_MATPLOTLIB: return
    types = ["Exon", "Intron+Exon", "Intron"]
    colors = ["#1A5084", "#914614", "#118730"]
    
    gene_data = [
        [stats_exon.get(w, {}).get("genes", 0) for w in wells],
        [stats_inex.get(w, {}).get("genes", 0) for w in wells],
        [stats_intron.get(w, {}).get("genes", 0) for w in wells],
    ]
    umi_data = [
        [stats_exon.get(w, {}).get("umis", 0) for w in wells],
        [stats_inex.get(w, {}).get("umis", 0) for w in wells],
        [stats_intron.get(w, {}).get("umis", 0) for w in wells],
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    
    bp0 = axes[0].boxplot(gene_data, notch=True, patch_artist=True)
    axes[0].set_xticklabels(types)
    for patch, color in zip(bp0['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axes[0].set_title("Number of Genes")
    axes[0].set_ylabel("Count")
    
    bp1 = axes[1].boxplot(umi_data, notch=True, patch_artist=True)
    axes[1].set_xticklabels(types)
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axes[1].set_title("Number of UMIs")
    axes[1].set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def calculate_saturation(out_dir, project):
    """
    Calculates saturation metrics:
    1. Seq_Saturation_Library (Global UMI distribution)
    2. Seq_Saturation_Gene (Gene-mapped UMI distribution)
    3. Median Genes per Cell (UMI based)
    4. Median Genes per Cell (Read based)
    """
    sat_json_global = os.path.join(stats_dir(out_dir), f"{project}.saturation_dist.json")
    sat_json_gene = os.path.join(stats_dir(out_dir), f"{project}.gene_saturation_dist.json")
    
    # Locate UMI and Read Matrices (Strictly enforce same type for comparison)
    umi_matrix_dir = None
    read_matrix_dir = None
    
    # Prefer Intron+Exon (inex) if available, else Exon
    for seq_type in ["inex", "exon"]:
        d_umi = os.path.join(expression_dir(out_dir), f"{project}.{seq_type}.umi")
        d_read = os.path.join(expression_dir(out_dir), f"{project}.{seq_type}.read")
        
        # Check for UMI matrix existence
        if os.path.exists(os.path.join(d_umi, "matrix.mtx.gz")):
            umi_matrix_dir = d_umi
            # Only use Read matrix if it matches the SAME type
            if os.path.exists(os.path.join(d_read, "matrix.mtx.gz")):
                read_matrix_dir = d_read
            else:
                read_matrix_dir = None # Reset to ensure we don't mix types
            
            print(f"Using sequence type '{seq_type}' for saturation/stats (UMI: Yes, Read: {'Yes' if read_matrix_dir else 'No'})")
            break
            
    if not os.path.exists(sat_json_global):
        print("Saturation distribution file not found. Skipping saturation analysis.")
        return None

    fractions = np.array([0.0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    def calc_sat_curve(json_file, fractions):
        if not os.path.exists(json_file): return np.zeros(len(fractions))
        
        with open(json_file, 'r') as f:
            umi_counts_dist = json.load(f)
        
        if isinstance(umi_counts_dist, dict):
            umi_counts = np.array(list(umi_counts_dist.keys()), dtype=np.int64)
            umi_freqs = np.array(list(umi_counts_dist.values()), dtype=np.float64)
        else:
            umi_counts = np.array(umi_counts_dist, dtype=np.int64)
            umi_freqs = np.ones_like(umi_counts, dtype=np.float64)

        total_umi_reads = np.sum(umi_counts * umi_freqs)
        if total_umi_reads == 0:
            return np.zeros(len(fractions))
        
        if umi_counts.size <= 2_000_000:
            probs = 1.0 - np.power(1.0 - fractions[:, np.newaxis], umi_counts[np.newaxis, :])
            exp_unique = np.dot(probs, umi_freqs)
            exp_reads = fractions * total_umi_reads
            saturation = np.where(exp_reads > 0, 1.0 - (exp_unique / exp_reads), 0.0)
            return saturation

        chunk_size = 100000
        n_chunks = (umi_counts.size + chunk_size - 1) // chunk_size
        exp_unique = np.zeros(len(fractions), dtype=np.float64)
        
        for i in range(n_chunks):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, umi_counts.size)
            chunk_counts = umi_counts[start:end]
            chunk_freqs = umi_freqs[start:end]
            probs = 1.0 - np.power(1.0 - fractions[:, np.newaxis], chunk_counts[np.newaxis, :])
            exp_unique += np.dot(probs, chunk_freqs)
        
        exp_reads = fractions * total_umi_reads
        saturation = np.where(exp_reads > 0, 1.0 - (exp_unique / exp_reads), 0.0)
        return saturation

    print(f"Calculating Library Saturation from {sat_json_global}...")
    seq_sat_lib = calc_sat_curve(sat_json_global, fractions)
    
    print(f"Calculating Gene Saturation from {sat_json_gene}...")
    seq_sat_gene = calc_sat_curve(sat_json_gene, fractions)

    # Helper to load matrix into per-cell counts
    def load_per_cell_counts(matrix_dir):
        if not matrix_dir: return None, 0
        m_path = os.path.join(matrix_dir, "matrix.mtx.gz")
        cell_data = collections.defaultdict(list)
        num_cells = 0
        try:
            opener = gzip.open if m_path.endswith('.gz') else open
            with opener(m_path, 'rt') as f:
                header_passed = False
                for line in f:
                    if line.startswith('%'): continue
                    if not header_passed:
                        parts = line.split()
                        num_cells = int(parts[1])
                        header_passed = True
                        continue
                    p = line.split()
                    cell_idx = int(p[1])
                    val = int(p[2])
                    if val > 0:
                        cell_data[cell_idx].append(val)
            return cell_data, num_cells
        except Exception as e:
            print(f"Error loading matrix {matrix_dir}: {e}")
            return None, 0

    def calc_median_gene_curve(cell_arrays, fractions):
        cell_arrays = [np.asarray(counts, dtype=np.float64) for counts in cell_arrays if len(counts)]
        if not cell_arrays:
            return [0.0] * len(fractions)

        lengths = np.fromiter((arr.size for arr in cell_arrays), dtype=np.int64)
        starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
        counts = np.concatenate(cell_arrays)

        medians = []
        for frac in fractions:
            probs = 1.0 - np.power(1.0 - frac, counts)
            per_cell = np.add.reduceat(probs, starts)
            medians.append(float(np.median(per_cell)))
        return medians

    # 2. Median Genes (UMI)
    median_genes_umi = []
    if umi_matrix_dir:
        print(f"Loading UMI Matrix for Median Genes: {umi_matrix_dir}")
        umi_cell_data, n_cells_u = load_per_cell_counts(umi_matrix_dir)
        print(f"  - Loaded {len(umi_cell_data)} cells (Header says {n_cells_u})")
        if umi_cell_data:
            print("Calculating Median Genes (UMI)...")
            # Optimize: Convert lists to numpy arrays once
            cell_arrays = [np.array(counts, dtype=np.int64) for counts in umi_cell_data.values()]
            
            # Debug stats
            total_counts = [c.sum() for c in cell_arrays]
            print(f"  - UMI Counts per cell: Mean={np.mean(total_counts):.1f}, Median={np.median(total_counts):.1f}")
            
            median_genes_umi = calc_median_gene_curve(cell_arrays, fractions)

    # 3. Median Genes (Read)
    median_genes_read = []
    if read_matrix_dir:
        print(f"Loading Read Matrix for Median Genes: {read_matrix_dir}")
        read_cell_data, n_cells_r = load_per_cell_counts(read_matrix_dir)
        print(f"  - Loaded {len(read_cell_data)} cells (Header says {n_cells_r})")
        if read_cell_data:
            print("Calculating Median Genes (Read)...")
            cell_arrays = [np.array(counts, dtype=np.int64) for counts in read_cell_data.values()]
            
            # Debug stats
            total_counts = [c.sum() for c in cell_arrays]
            print(f"  - Read Counts per cell: Mean={np.mean(total_counts):.1f}, Median={np.median(total_counts):.1f}")

            median_genes_read = calc_median_gene_curve(cell_arrays, fractions)
            
    return {
        "fractions": fractions,
        "seq_saturation_lib": seq_sat_lib,
        "seq_saturation_gene": seq_sat_gene,
        "median_genes_umi": median_genes_umi,
        "median_genes_read": median_genes_read
    }

def plot_saturation(stats, out_dir, project):
    if not HAS_MATPLOTLIB or not stats: return
    
    fractions = stats['fractions']
    seq_sat_lib = stats['seq_saturation_lib']
    seq_sat_gene = stats['seq_saturation_gene']
    med_umi = stats.get('median_genes_umi')
    med_read = stats.get('median_genes_read')
    
    fig, ax1 = plt.subplots(figsize=(8, 6))
    
    # Plot Sequencing Saturation
    ax1.set_xlabel('Sequencing Depth (Fraction)')
    ax1.set_ylabel('Sequencing Saturation', color='black')
    
    ax1.plot(fractions, seq_sat_lib, marker='o', color='tab:red', label='Sat (Library)')
    ax1.plot(fractions, seq_sat_gene, marker='x', linestyle='--', color='tab:orange', label='Sat (Gene)')
    
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.set_ylim(0, 1.05)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Plot Median Genes
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Median Genes per Cell', color=color)
    
    if med_umi:
        ax2.plot(fractions, med_umi, marker='s', linestyle='-', color=color, label='Median Genes (UMI)')
    if med_read:
        ax2.plot(fractions, med_read, marker='^', linestyle='--', color='tab:green', label='Median Genes (Read)')
        
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Combined Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')
    
    plt.title(f"Saturation Analysis: {project}")
    fig.tight_layout()
    
    out_pdf = os.path.join(stats_dir(out_dir), f"{project}.saturation.pdf")
    plt.savefig(out_pdf)
    plt.close()
    
    # Write TSV
    out_tsv = os.path.join(stats_dir(out_dir), f"{project}.saturation.tsv")
    print(f"Writing saturation table to {out_tsv}...")
    with open(out_tsv, 'w') as f:
        # Updated Header
        f.write("Fraction\tSeq_Saturation_Library\tSeq_Saturation_Gene\tMedian_Genes_UMI\tMedian_Genes_Read\n")
        
        for i, frac in enumerate(fractions):
            val_umi = med_umi[i] if med_umi and i < len(med_umi) else 0.0
            val_read = med_read[i] if med_read and i < len(med_read) else 0.0
            
            f.write(f"{frac:.1f}\t{seq_sat_lib[i]:.4f}\t{seq_sat_gene[i]:.4f}\t{val_umi:.2f}\t{val_read:.2f}\n")

def plot_features(read_stats, wells, out_pdf):
    if not HAS_MATPLOTLIB: return
    feat_colors = {
        "Exon": "#1A5084", "Intron": "#118730", "Unmapped": "#545454",
        "Ambiguity": "#FFA54F", 
        "Intergenic": "#FFD700", "Other_Unassigned": "#8A8F98", "Unused BC": "#BABABA",
    }
    categories_order = ["Exon", "Intron", "Intergenic", "Ambiguity", "Other_Unassigned", "Unmapped"]
    
    totals = {c: 0 for c in categories_order}
    per_cell_frac = {c: [] for c in categories_order}
    unused_total = int(read_stats.get("__NO_CB__", {}).get("Unused BC", 0))

    for bc in wells:
        counts = read_stats.get(bc, {})
        cell_total = 0
        for c in categories_order:
            v = int(counts.get(c, 0))
            totals[c] += v
            cell_total += v
        if cell_total > 0:
            for c in categories_order:
                per_cell_frac[c].append(int(counts.get(c, 0)) / cell_total)
        else:
            for c in categories_order: per_cell_frac[c].append(0.0)

    total_sum = sum(totals.values()) + unused_total
    total_fracs = {c: (totals[c] / total_sum if total_sum > 0 else 0.0) for c in categories_order}
    unused_frac = unused_total / total_sum if total_sum > 0 else 0.0

    bar_categories = [c for c in categories_order if totals.get(c, 0) > 0]
    box_categories = [c for c in categories_order if len(per_cell_frac[c]) > 0]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [0.6, 1.0]})

    left = 0.0
    for c in bar_categories:
        w = total_fracs.get(c, 0.0)
        if w > 0:
            axes[0].barh([0], [w * 100.0], left=[left * 100.0], color=feat_colors[c], label=c, height=0.5)
            left += w
    if unused_frac > 0:
        axes[0].barh([0], [unused_frac * 100.0], left=[left * 100.0], color=feat_colors["Unused BC"], label="Unused BC", height=0.5)

    axes[0].set_xlim(0, 100)
    axes[0].set_yticks([])
    axes[0].set_xlabel("% of total reads")
    axes[0].set_title("Total Read Distribution")
    axes[0].legend(ncol=max(1, len(bar_categories)+1), loc="upper center", bbox_to_anchor=(0.5, -0.30))

    box_pairs = []
    for c in box_categories:
        vals = [v * 100.0 for v in per_cell_frac[c]]
        if len(vals) == 0:
            continue
        box_pairs.append((c, vals))

    if not box_pairs or total_sum <= 0:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No reads available for per-cell distribution", ha="center", va="center")
    else:
        plot_labels = [c for c, _vals in box_pairs]
        box_data = [_vals for _c, _vals in box_pairs]
        bp = axes[1].boxplot(box_data, notch=True, patch_artist=True)
        axes[1].set_xticklabels(plot_labels)
        for patch, c in zip(bp["boxes"], plot_labels):
            patch.set_facecolor(feat_colors[c])

    axes[1].set_title("Reads per Cell")
    axes[1].set_ylabel("% reads/cell")
    
    fig.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 generate_stats.py <yaml_config>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    
    stats_output_dir = stats_dir(out_dir)
    if not os.path.exists(stats_output_dir): os.makedirs(stats_output_dir)
    
    # 1. Load Barcode Mapping
    bc_mapping = load_barcode_mapping(out_dir, project)
    all_wells = sorted(list(bc_mapping.keys()))
    kept_barcodes = set(all_wells)
    
    # 2. Calculate Matrix Stats (Fast, independent)
    print("Calculating Matrix Statistics...")
    umi_exon_dir = os.path.join(expression_dir(out_dir), f"{project}.exon.umi")
    umi_intron_dir = os.path.join(expression_dir(out_dir), f"{project}.intron.umi")
    
    compact_stats = load_cell_matrix_stats(os.path.join(stats_output_dir, f"{project}.cell_matrix_stats.json"))

    # Load Read stats (New Feature)
    read_exon_dir = os.path.join(expression_dir(out_dir), f"{project}.exon.read")
    read_intron_dir = os.path.join(expression_dir(out_dir), f"{project}.intron.read")

    if compact_stats is not None:
        print("Using DGE cell matrix QC summary; skipping repeated MEX parsing.")
        stats_exon = compact_stats["umi"].get("exon", {})
        stats_intron = compact_stats["umi"].get("intron", {})
        stats_inex = compact_stats["umi"].get("inex", {})
        rstats_exon = compact_stats["read"].get("exon", {})
        rstats_intron = compact_stats["read"].get("intron", {})
        rstats_inex = compact_stats["read"].get("inex", {})
    else:
        stats_exon = calculate_matrix_stats(umi_exon_dir)
        stats_intron = calculate_matrix_stats(umi_intron_dir)
        stats_inex = collections.defaultdict(lambda: {'umis': 0, 'genes': 0})
        for bc in set(stats_exon.keys()) | set(stats_intron.keys()):
            ex_st = stats_exon.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
            in_st = stats_intron.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
            stats_inex[bc]['umis'] = ex_st['umis'] + in_st['umis']
            stats_inex[bc]['genes'] = len(ex_st['gene_set'] | in_st['gene_set'])

        rstats_exon = calculate_matrix_stats(read_exon_dir)
        rstats_intron = calculate_matrix_stats(read_intron_dir)
        rstats_inex = collections.defaultdict(lambda: {'genes': 0})
        for bc in set(rstats_exon.keys()) | set(rstats_intron.keys()):
            ex_st = rstats_exon.get(bc, {'genes': 0, 'gene_set': set()})
            in_st = rstats_intron.get(bc, {'genes': 0, 'gene_set': set()})
            rstats_inex[bc]['genes'] = len(ex_st['gene_set'] | in_st['gene_set'])

    # 3. Load Pre-calculated Read Stats & Coverage
    stats_json_path = os.path.join(stats_output_dir, f"{project}.read_stats.json")
    
    read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_umi = np.zeros(100, dtype=np.float64)
    cov_int = np.zeros(100, dtype=np.float64)
    
    if os.path.exists(stats_json_path):
        print(f"Loading pre-calculated stats from {stats_json_path}")
        with open(stats_json_path, 'r') as f:
            data = json.load(f)
            # Reconstruct defaultdict from JSON dict
            raw_stats = data.get("read_stats", {})
            for k, v in raw_stats.items():
                read_stats[k] = v
                
            cov_umi_list = data.get("coverage_umi", [])
            cov_int_list = data.get("coverage_int", [])
            
            if cov_umi_list:
                cov_umi = np.array(cov_umi_list, dtype=np.float64)
            if cov_int_list:
                cov_int = np.array(cov_int_list, dtype=np.float64)
    else:
        print("Warning: Pre-calculated stats not found. Plots will be missing.")

    # 4. Plots and Tables
    out_prefix = os.path.join(stats_output_dir, project)
    plot_coverage(cov_umi, cov_int, out_prefix)
    
    # Write Table
    output_table = os.path.join(stats_output_dir, f"{project}.stats.tsv")
    print(f"Writing stats table to {output_table}...")
    
    with open(output_table, 'w') as f:
        headers = [
            "wellID", "internal_barcodes", "umi_barcodes", 
            "internal_reads", "umi_reads", "all_reads",
            "Ambiguity_reads", "Exon_reads", "Intergenic_reads", "intron_reads", 
            "Unmapped_reads", "Other_Unassigned_reads", "MappingRatio", "GenicRatio", "ExonIntronRatio", "UMIfrac",
            "Exon_umis", "Intron_umis", "Intron_Exon_umis",
            "Exon_genes", "Intron_genes", "Intron_Exon_genes",
            "Exon_read_genes", "Intron_read_genes", "Intron_Exon_read_genes"
        ]
        f.write("\t".join(headers) + "\n")
        
        for well in all_wells:
            int_bc = bc_mapping[well]['int']
            umi_bc = bc_mapping[well]['umi']
            
            r_stats = read_stats.get(well, {})
            ambig = r_stats.get('Ambiguity', 0)
            exon_r = r_stats.get('Exon', 0)
            inter_r = r_stats.get('Intergenic', 0)
            intron_r = r_stats.get('Intron', 0)
            unmapped = r_stats.get('Unmapped', 0)
            other_unassigned = r_stats.get('Other_Unassigned', 0)

            internal_r = r_stats.get('Internal_Reads', 0)
            umi_r = r_stats.get('UMI_Reads', 0)
            all_r = internal_r + umi_r
            mapping_ratio, genic_ratio, exon_intron_ratio = calculate_read_ratios(
                exon_r,
                intron_r,
                inter_r,
                ambig,
                unmapped,
                other_unassigned,
                all_r,
            )
            umi_frac = (umi_r / all_r) if all_r > 0 else ""
            
            ex_st = stats_exon.get(well, {'umis': 0, 'genes': 0})
            in_st = stats_intron.get(well, {'umis': 0, 'genes': 0})
            inex_st = stats_inex.get(well, {'umis': 0, 'genes': 0})
            
            rex_st = rstats_exon.get(well, {'genes': 0})
            rin_st = rstats_intron.get(well, {'genes': 0})
            rinex_st = rstats_inex.get(well, {'genes': 0})
            
            row = [
                well, int_bc, umi_bc,
                internal_r, umi_r, all_r,
                ambig, exon_r, inter_r, intron_r,
                unmapped, other_unassigned,
                f"{mapping_ratio:.6f}" if mapping_ratio is not None else "",
                f"{genic_ratio:.6f}" if genic_ratio is not None else "",
                f"{exon_intron_ratio:.6f}" if exon_intron_ratio is not None else "",
                f"{umi_frac:.6f}" if umi_frac != "" else "",
                ex_st['umis'], in_st['umis'], inex_st['umis'],
                ex_st['genes'], in_st['genes'], inex_st['genes'],
                rex_st['genes'], rin_st['genes'], rinex_st['genes']
            ]
            f.write("\t".join(map(str, row)) + "\n")

    if HAS_MATPLOTLIB:
        plot_gene_umi_counts_by_type(stats_exon, stats_intron, stats_inex, all_wells, os.path.join(stats_output_dir, f"{project}.geneUMIcounts.pdf"))
        plot_features(read_stats, kept_barcodes, os.path.join(stats_output_dir, f"{project}.features.pdf"))
            
    # 5. Saturation Analysis
    sat_stats = calculate_saturation(out_dir, project)
    if sat_stats:
        plot_saturation(sat_stats, out_dir, project)
        
    print("Statistics generation finished.")

if __name__ == "__main__":
    main()
