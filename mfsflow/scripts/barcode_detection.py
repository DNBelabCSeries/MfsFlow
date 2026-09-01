#!/usr/bin/env python3
"""
Barcode detection and processing for standalone script execution.

This module handles barcode detection from sequencing data, including
whitelist loading, barcode counting, and statistical analysis for
quality assessment of single-cell experiments.
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys
import re
import collections
import multiprocessing as mp

try:
    from mfsflow.path_layout import barcode_dir, load_config
except ImportError:
    from path_layout import barcode_dir, load_config

# Suppress pandas chained assignment warnings
pd.options.mode.chained_assignment = None

def read_whitelist(bcfile):
    """
    Reads a whitelist file robustly.
    Handles:
    - Standard 1-column or 2-column files.
    - Comma-separated barcodes within columns/rows.
    - Mixed whitespace delimiters (tabs, spaces).
    """
    if not os.path.exists(bcfile):
        print(f"Error: Whitelist file {bcfile} not found.")
        sys.exit(1)
        
    try:
        with open(bcfile, 'r') as f:
            content = f.read()
        
        # Robust splitting: treats commas, tabs, spaces, newlines all as delimiters
        tokens = re.split(r'[,\s]+', content)
        # Remove empty strings resulting from consecutive delimiters
        bc_wl = sorted(list(set([t for t in tokens if t])))
        
        return set(bc_wl)
    except Exception as e:
        print(f"Error reading whitelist: {e}")
        sys.exit(1)


def write_unmatched_whitelist_barcodes(bccount_df, whitelist, out_file, limit=100):
    """Write the most abundant observed barcodes absent from the whitelist."""
    limit = max(0, int(limit))
    unmatched = bccount_df.loc[~bccount_df['XC'].isin(whitelist), ['XC', 'n']].copy()
    unmatched = unmatched.sort_values(by=['n', 'XC'], ascending=[False, True]).head(limit)
    unmatched = unmatched.rename(columns={'XC': 'barcode', 'n': 'reads'})

    parent = os.path.dirname(out_file)
    if parent:
        os.makedirs(parent, exist_ok=True)
    unmatched.to_csv(out_file, sep='\t', index=False)
    return len(unmatched)


def load_expected_id_map(expect_file):
    expected = collections.OrderedDict()
    if not expect_file or not os.path.exists(expect_file):
        return expected

    with open(expect_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0].lower() == "wellid":
                continue
            well_id = parts[0]
            umi_barcodes = [x.strip().upper() for x in parts[1].split(",") if x.strip()]
            internal_barcodes = [x.strip().upper() for x in parts[2].split(",") if x.strip()]
            expected[well_id] = {
                "umi": umi_barcodes,
                "internal": internal_barcodes,
            }
    return expected


def write_barcodes_by_well(kept_df, expect_file, out_file):
    expected = load_expected_id_map(expect_file)
    if not expected:
        return False

    count_map = {
        str(row.XC).strip().upper(): int(row.n)
        for row in kept_df.itertuples(index=False)
    }
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, "w") as out:
        out.write("wellID\tumi_reads\tinternal_reads\ttotal_reads\tumi_barcodes\tinternal_barcodes\n")
        for well_id, values in expected.items():
            umi_reads = sum(count_map.get(bc, 0) for bc in values["umi"])
            internal_reads = sum(count_map.get(bc, 0) for bc in values["internal"])
            total_reads = umi_reads + internal_reads
            out.write(
                f"{well_id}\t{umi_reads}\t{internal_reads}\t{total_reads}\t"
                f"{','.join(values['umi'])}\t{','.join(values['internal'])}\n"
            )
    return True

def find_knee_point(sorted_counts):
    """
    Geometric implementation of the Knee/Elbow method.
    Approximates the 'uik' method in R.
    """
    if len(sorted_counts) == 0:
        return 0
        
    y = np.log10(sorted_counts)
    x = np.arange(len(y))
    
    start_point = np.array([x[0], y[0]])
    end_point = np.array([x[-1], y[-1]])
    vec_line = end_point - start_point
    
    vec_points = np.column_stack((x, y)) - start_point
    
    cross_prod = vec_points[:, 0] * vec_line[1] - vec_points[:, 1] * vec_line[0]
    distances = np.abs(cross_prod) / np.linalg.norm(vec_line)
    
    knee_idx = np.argmax(distances)
    return sorted_counts[knee_idx]

def cell_bc_selection(bccount_df, config):
    """
    Select true cell barcodes using the configured whitelist/automatic strategy.
    Known whitelist modes are intentionally strict and do not fall back to top barcodes.
    """
    barcodes_config = config['barcodes']
    # Keep in sync with DEFAULT_CONFIG["barcodes"]["nReadsperCell"] (config/defaults.py).
    min_reads = barcodes_config.get('nReadsperCell', 1)
    
    # Keep the full barcode table so downstream Hamming rescue can still see
    # low-count tail barcodes. Selection itself is still thresholded.
    df = bccount_df.copy()
    df = df.sort_values(by=['n', 'XC'], ascending=[False, True]).reset_index(drop=True)
    df['keep'] = False
    eligible = df[df['n'] >= min_reads].copy()
    
    strategy_auto = barcodes_config.get('automatic', False)
    bc_num = barcodes_config.get('barcode_num', None)
    bc_file = barcodes_config.get('barcode_file', None)
    
    whitelist = set()
    if bc_file:
        whitelist = read_whitelist(bc_file)
        if not whitelist:
            raise ValueError(f"Whitelist file is empty: {bc_file}")

    # --- Strategy Selection ---
    
    if bc_file and strategy_auto:
        # Strategy: Automatic Knee + Whitelist Intersection
        print("Strategy: Automatic + Whitelist Intersection")
        cutoff = find_knee_point(eligible['n'].values)
        print(f"  Automatic cutoff: {cutoff} reads")
        
        # Mark potential cells by cutoff
        potential_keep = eligible['n'] >= cutoff
        
        # Filter those by whitelist
        # We strictly check if the potential high-read BC is in whitelist
        valid_in_whitelist = eligible.loc[potential_keep, 'XC'].isin(whitelist)
        
        if valid_in_whitelist.any():
            # If we have intersection, keep them
            selected = set(eligible.loc[potential_keep & eligible['XC'].isin(whitelist), 'XC'])
            df.loc[df['XC'].isin(selected), 'keep'] = True
        else:
            raise ValueError(
                "Automatic detection found no overlap with whitelist. "
                "Check barcode position or --manual/--plate selection."
            )

    elif bc_file and not strategy_auto:
        # Strategy: Strict Whitelist
        print("Strategy: Known Whitelist")
        observed_in_whitelist = df['XC'].isin(whitelist)
        eligible_in_whitelist = eligible['XC'].isin(whitelist)
        df['keep'] = df['XC'].isin(set(eligible.loc[eligible_in_whitelist, 'XC']))

        if not observed_in_whitelist.any():
            raise ValueError(
                "None of the annotated barcodes were detected. "
                "The pipeline will not fall back to top barcodes for manual/plate/custom modes; "
                "check barcode position or --manual/--plate selection."
            )

    elif bc_num is not None:
        # Strategy: Fixed Number
        print(f"Strategy: Fixed Number ({bc_num})")
        limit = min(int(bc_num), len(eligible))
        selected = set(eligible.iloc[:limit]['XC'])
        df.loc[df['XC'].isin(selected), 'keep'] = True
        
    else:
        # Strategy: Automatic (Knee)
        print("Strategy: Automatic (Knee method)")
        cutoff = find_knee_point(eligible['n'].values)
        print(f"  Automatic cutoff: {cutoff} reads")
        selected = set(eligible.loc[eligible['n'] >= cutoff, 'XC'])
        df.loc[df['XC'].isin(selected), 'keep'] = True
        
        # Fallback for too few cells
        if df['keep'].sum() < 10:
             print("  Warning: < 10 cells found. Using top 100 fallback.")
             limit = min(100, len(eligible))
             selected = set(eligible.iloc[:limit]['XC'])
             df.loc[df['XC'].isin(selected), 'keep'] = True

    print(f"Selected {df['keep'].sum()} cell barcodes.")
    return df

_HB_TRUE_SET = None
_HB_MASK_TO_TRUE = None
_HB_BC_LEN = None

def _hb_init(true_set, mask_to_true, bc_len):
    global _HB_TRUE_SET, _HB_MASK_TO_TRUE, _HB_BC_LEN
    _HB_TRUE_SET = true_set
    _HB_MASK_TO_TRUE = mask_to_true
    _HB_BC_LEN = bc_len

def _hb_process_chunk(cand_chunk):
    out = []
    true_set = _HB_TRUE_SET
    mask_to_true = _HB_MASK_TO_TRUE
    bc_len = _HB_BC_LEN

    for cand in cand_chunk:
        if len(cand) != bc_len:
            continue
        if cand in true_set:
            out.append((cand, cand, 0))
            continue
        matches = set()
        for i in range(bc_len):
            mask = cand[:i] + '*' + cand[i+1:]
            for t in mask_to_true.get(mask, ()):
                matches.add(t)
        if matches:
            for t in matches:
                out.append((cand, t, 1))
    return out

def fast_hamming_binning(true_bcs, candidate_bcs, threshold=1, threads=1):
    """
    Optimized Hamming distance binning.
    Returns: (final_df, raw_df)
    - raw_df: All matches <= threshold (ties included)
    - final_df: Unambiguous matches only (ties discarded)
    """
    true_list = [str(x).strip().upper() for x in true_bcs if pd.notna(x)]
    cand_list = [str(x).strip().upper() for x in candidate_bcs if pd.notna(x)]

    if not true_list or not cand_list:
        return pd.DataFrame(), pd.DataFrame()

    true_set = set(true_list)
    bc_len = len(next(iter(true_set)))

    if threshold < 0:
        return pd.DataFrame(), pd.DataFrame()

    mapping_columns = ["falseBC", "trueBC", "hamming"]
    mapping_list = []

    if threshold == 0:
        for cand in cand_list:
            if len(cand) != bc_len:
                continue
            if cand in true_set:
                mapping_list.append((cand, cand, 0))
        raw_df = pd.DataFrame(mapping_list, columns=mapping_columns)
    elif threshold == 1:
        mask_to_true = collections.defaultdict(list)
        for t in true_set:
            if len(t) != bc_len:
                continue
            for i in range(bc_len):
                mask_to_true[t[:i] + '*' + t[i+1:]].append(t)

        threads = int(threads) if threads else 1
        if threads > 1 and len(cand_list) >= 20000:
            print(f"Starting binning: {len(true_set)} True BCs vs {len(cand_list)} Candidate BCs (threshold=1, threads={threads})")
            try:
                ctx = mp.get_context('fork')
            except ValueError:
                ctx = mp.get_context()

            chunk_size = 5000
            cand_chunks = [cand_list[i:i+chunk_size] for i in range(0, len(cand_list), chunk_size)]

            with ctx.Pool(
                processes=threads,
                initializer=_hb_init,
                initargs=(true_set, dict(mask_to_true), bc_len),
            ) as pool:
                for i, out in enumerate(pool.imap_unordered(_hb_process_chunk, cand_chunks), start=1):
                    mapping_list.extend(out)
                    if i % 10 == 0 or i == len(cand_chunks):
                        print(f"  Processed chunks: {i}/{len(cand_chunks)}", end='\r')
            print()
        else:
            print(f"Starting binning: {len(true_set)} True BCs vs {len(cand_list)} Candidate BCs (threshold=1)")
            for cand in cand_list:
                if len(cand) != bc_len:
                    continue
                if cand in true_set:
                    mapping_list.append((cand, cand, 0))
                    continue
                matches = set()
                for i in range(bc_len):
                    mask = cand[:i] + '*' + cand[i+1:]
                    for t in mask_to_true.get(mask, ()):
                        matches.add(t)
                if matches:
                    for t in matches:
                        mapping_list.append((cand, t, 1))
        raw_df = pd.DataFrame(mapping_list, columns=mapping_columns)
    else:
        print(f"Warning: BarcodeBinning={threshold} may be slow. Consider using 0 or 1 for performance.")
        true_arr = np.array(sorted(true_set), dtype=object)
        chunk_size = 2000
        total_chunks = (len(cand_list) // chunk_size) + 1
        print(f"Starting binning: {len(true_arr)} True BCs vs {len(cand_list)} Candidate BCs (threshold={threshold})")
        for i in range(0, len(cand_list), chunk_size):
            if i % (chunk_size * 10) == 0:
                print(f"  Processing chunk {i // chunk_size + 1}/{total_chunks}...")
            for cand in cand_list[i:i+chunk_size]:
                if len(cand) != bc_len:
                    continue
                best = threshold + 1
                best_matches = []
                for t in true_arr:
                    dist = 0
                    for c1, c2 in zip(cand, t):
                        if c1 != c2:
                            dist += 1
                            if dist > best:
                                break
                    if dist < best:
                        best = dist
                        best_matches = [t]
                    elif dist == best:
                        best_matches.append(t)
                if best <= threshold:
                    for t in best_matches:
                        mapping_list.append((cand, t, best))
        raw_df = pd.DataFrame(mapping_list, columns=mapping_columns)
    
    if raw_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Filter for Ambiguity (Matches R logic)
    # R: binmap_raw[n_min==1]
    
    # 1. Count how many trueBCs each falseBC maps to (at min distance)
    # Since we only stored min_dist matches in mapping_list, we just count occurrences
    counts = raw_df.groupby('falseBC')['trueBC'].count()
    
    # 2. Keep only those with count == 1 (Unambiguous)
    valid_false_bcs = counts[counts == 1].index
    final_df = raw_df[raw_df['falseBC'].isin(valid_false_bcs)].copy()
    
    return final_df, raw_df

def main():
    parser = argparse.ArgumentParser(description="Barcode detection for MfsFlow")
    parser.add_argument('yaml_config', help="Path to run config YAML file")
    args = parser.parse_args()
    
    print(f"Loading config from {args.yaml_config}")
    opt = load_config(args.yaml_config)
    
    project = opt['project']
    out_dir = opt['out_dir']
    
    # Ensure output directory
    output_base = barcode_dir(out_dir)
    if not os.path.exists(output_base):
        os.makedirs(output_base)
        
    bccount_file = os.path.join(out_dir, f"{project}.BCstats.txt")
    print(f"Reading barcode stats from {bccount_file}...")
    
    if not os.path.exists(bccount_file):
        print(f"Error: Stats file {bccount_file} not found.")
        sys.exit(1)

    # Read raw counts
    raw_df = pd.read_csv(bccount_file, sep='\t', header=None, names=['XC', 'n'])
    # Aggregate counts just in case of duplicates
    raw_df = raw_df.groupby('XC', as_index=False)['n'].sum()

    # Write this diagnostic before selection so it remains available when a
    # strict whitelist has no overlap and cell_bc_selection raises an error.
    unmatched_file = None
    bc_file = opt.get('barcodes', {}).get('barcode_file')
    if bc_file:
        whitelist = read_whitelist(bc_file)
        unmatched_file = os.path.join(output_base, f"{project}.unmatched_whitelist_barcodes.tsv")
        unmatched_count = write_unmatched_whitelist_barcodes(
            raw_df,
            whitelist,
            unmatched_file,
            limit=100,
        )
        print(
            f"Wrote {unmatched_count} unmatched whitelist barcode(s) "
            f"(top 100 by reads): {unmatched_file}"
        )
    
    # --- Step 1: Selection ---
    try:
        df_processed = cell_bc_selection(raw_df, opt)
    except ValueError as exc:
        if unmatched_file:
            raise ValueError(f"{exc} Top unmatched barcodes: {unmatched_file}") from exc
        raise
    
    # Filter to kept only
    kept_df = df_processed[df_processed['keep']].copy()
    
    # Save Step 1 Result (Pre-binning)
    # Rfwrite defaults to including headers "XC" and "n"
    kept_file = os.path.join(output_base, f"{project}kept_barcodes.txt")
    kept_df[['XC', 'n']].to_csv(kept_file, sep='\t', index=False)
    empty_bin_file = os.path.join(output_base, f"{project}.BCbinning.txt")
    
    # --- Step 2: Binning ---
    bc_opts = opt.get('barcodes', {})
    do_binning = False
    
    if bc_opts.get('BarcodeBinning', 0) > 0:
        do_binning = True
    if bc_opts.get('barcode_sharing'):
        # Generic barcode sharing uses the same Hamming recovery path here.
        # The original substring replacement table logic is not enabled by default.
        do_binning = True
        print("Note: generic barcode sharing requested; using Hamming recovery path.")

    if do_binning and len(kept_df) > 0:
        true_bcs = kept_df['XC'].values
        # Candidates: all non-kept barcodes, including low-count tail barcodes.
        candidates_df = df_processed[~df_processed['keep']]
        candidates_bcs = candidates_df['XC'].values
        
        if len(candidates_bcs) > 0:
            threshold = bc_opts.get('BarcodeBinning', 1)
            bin_map, bin_map_raw = fast_hamming_binning(true_bcs, candidates_bcs, threshold=threshold, threads=opt.get('num_threads', 1))
            
            if not bin_map.empty:
                print(f"Binned {len(bin_map)} barcodes.")
                
                # Pre-calculate counts map
                candidates_df_indexed = candidates_df.set_index('XC')
                
                # --- Save Raw Map ---
                raw_file = os.path.join(output_base, f"{project}.BCbinning.raw.txt")
                bin_map_raw['n'] = bin_map_raw['falseBC'].map(candidates_df_indexed['n']).fillna(0).astype(int)
                # Raw file usually needs specific columns too? R just dumps data.table
                # We align with binmap structure: falseBC, hamming, trueBC, n
                bin_map_raw = bin_map_raw[['falseBC', 'hamming', 'trueBC', 'n']]
                bin_map_raw.to_csv(raw_file, sep=',', index=False)

                # --- Save Final Map ---
                bin_file = os.path.join(output_base, f"{project}.BCbinning.txt")
                
                # Get 'n' for falseBCs (already done for raw, reuse logic if needed, or just copy)
                # bin_map is a subset of raw, so we can just filter raw or re-map
                bin_map['n'] = bin_map['falseBC'].map(candidates_df_indexed['n']).fillna(0).astype(int)
                
                bin_map_output = bin_map[['falseBC', 'hamming', 'trueBC', 'n']]
                bin_map_output.to_csv(bin_file, sep=',', index=False)
                
                # --- Update Counts ---
                add_counts = bin_map.groupby('trueBC')['n'].sum()
                
                kept_df = kept_df.set_index('XC')
                kept_df['n'] = kept_df['n'].add(add_counts, fill_value=0)
                kept_df = kept_df.reset_index()
                
                # Save Final Result
                binned_file = os.path.join(output_base, f"{project}kept_barcodes_binned.txt")
                kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
            else:
                print("No barcodes binned (no matches within threshold).")
                pd.DataFrame(columns=['falseBC', 'hamming', 'trueBC', 'n']).to_csv(empty_bin_file, sep=',', index=False)
                binned_file = os.path.join(output_base, f"{project}kept_barcodes_binned.txt")
                kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
        else:
            print("No candidate barcodes for binning.")
            pd.DataFrame(columns=['falseBC', 'hamming', 'trueBC', 'n']).to_csv(empty_bin_file, sep=',', index=False)
            # Save copy
            binned_file = os.path.join(output_base, f"{project}kept_barcodes_binned.txt")
            kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
    
    elif do_binning and len(kept_df) == 0:
        print("Warning: No true barcodes selected. Skipping binning.")
    
    print("Done.")

    expect_file = os.path.join(out_dir, "config", "expect_id_barcode.tsv")
    by_well_file = os.path.join(output_base, f"{project}.kept_barcodes_by_well.tsv")
    if write_barcodes_by_well(kept_df, expect_file, by_well_file):
        print(f"Wrote barcode reads by well: {by_well_file}")

if __name__ == "__main__":
    main()
