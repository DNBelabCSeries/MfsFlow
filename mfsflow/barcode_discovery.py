"""
Barcode discovery and expected barcode table generation for sample identification.

This module provides functions to build expected barcode records from manual
and auto barcode lists, write barcode tables for the pipeline, and perform
barcode discovery from sequencing data statistics.
"""

import os
import collections

import yaml


def build_expected_records(script_dir, sample_type=None, sample_ids=None):
    """Build expected barcode records from manual and/or auto barcode lists.
    
    Args:
        script_dir (str): Directory containing barcode YAML files.
        sample_type (str, optional): Sample type filter ('manual', 'auto', or None for both).
        sample_ids (str or list, optional): Sample IDs to filter.
        
    Returns:
        list: List of barcode record dictionaries.
    """
    records = []
    if sample_type in (None, "manual"):
        records.extend(_manual_records(script_dir, sample_ids if sample_type == "manual" else None))
    if sample_type in (None, "auto"):
        records.extend(_auto_records(script_dir, sample_ids if sample_type == "auto" else None))
    return records


def write_expected_tables(records, config_dir):
    """Write expected barcode tables for pipeline consumption.
    
    Args:
        records (list): List of barcode record dictionaries.
        config_dir (str): Directory to write output files.
        
    Returns:
        tuple: Paths to (pipe_path, summary_path).
        
    Raises:
        ValueError: If no records are provided.
    """
    if not records:
        raise ValueError("No expected barcodes were generated.")

    grouped = collections.OrderedDict()
    for rec in records:
        well_id = rec["wellID"]
        if well_id not in grouped:
            grouped[well_id] = {"umi": [], "internal": []}
        grouped[well_id][rec["barcode_type"]].append(rec["barcode"])

    summary_path = os.path.join(config_dir, "expect_id_barcode.tsv")
    pipe_path = os.path.join(config_dir, "expect_barcode.tsv")
    with open(summary_path, "w") as summary, open(pipe_path, "w") as pipe:
        print("wellID\tumi_barcodes\tinternal_barcodes", file=summary)
        for well_id, values in grouped.items():
            umi = _unique(values["umi"])
            internal = _unique(values["internal"])
            print(f"{well_id}\t{','.join(umi)}\t{','.join(internal)}", file=summary)
            for barcode in umi + internal:
                print(barcode, file=pipe)

    return pipe_path, summary_path


def discover_barcodes(bcstats_file, records, out_file, max_hamming=1, min_unique_barcodes=2, min_fraction=0.2):
    """Discover barcode candidates from observed barcode statistics.
    
    Compares observed barcodes against expected records, allowing for
    Hamming distance mismatches, and selects the most confident candidates.
    
    Args:
        bcstats_file (str): Path to barcode statistics file.
        records (list): List of expected barcode record dictionaries.
        out_file (str): Path to write discovery report.
        max_hamming (int, optional): Maximum Hamming distance for matching. Defaults to 1.
        min_unique_barcodes (int, optional): Minimum unique barcodes for a candidate. Defaults to 2.
        min_fraction (float, optional): Minimum fraction of reads for a candidate. Defaults to 0.2.
        
    Returns:
        tuple: (selected_candidates, selected_records).
        
    Raises:
        ValueError: If no barcode counts found or no confident candidates.
    """
    observed = _read_bcstats(bcstats_file)
    if not observed:
        raise ValueError(f"No barcode counts found in {bcstats_file}")
    return discover_barcodes_from_counts(
        observed,
        records,
        out_file,
        max_hamming=max_hamming,
        min_unique_barcodes=min_unique_barcodes,
        min_fraction=min_fraction,
    )


def discover_barcodes_from_counts(
    observed,
    records,
    out_file,
    max_hamming=1,
    min_unique_barcodes=2,
    min_fraction=0.2,
):
    """Discover barcode candidates from observed barcode counts."""
    if not observed:
        raise ValueError("No barcode counts found")

    exact = collections.defaultdict(list)
    masked = collections.defaultdict(list)
    for rec in records:
        barcode = rec["barcode"]
        exact[barcode].append(rec)
        if max_hamming >= 1:
            for i in range(len(barcode)):
                masked[barcode[:i] + "*" + barcode[i + 1:]].append(rec)

    candidate_stats = collections.defaultdict(lambda: {
        "candidate_type": "",
        "candidate_id": "",
        "matched_reads": 0,
        "exact_reads": 0,
        "hamming1_reads": 0,
        "matched_observed_barcodes": set(),
        "matched_expected_barcodes": set(),
    })
    match_rows = []

    for obs_bc, reads in observed.items():
        matches = []
        if obs_bc in exact:
            matches = [(0, rec) for rec in exact[obs_bc]]
        elif max_hamming >= 1:
            seen = set()
            for i in range(len(obs_bc)):
                key = obs_bc[:i] + "*" + obs_bc[i + 1:]
                for rec in masked.get(key, ()):
                    rec_key = (rec["candidate_type"], rec["candidate_id"], rec["wellID"], rec["barcode"], rec["barcode_type"])
                    if rec_key not in seen:
                        seen.add(rec_key)
                        matches.append((1, rec))

        if not matches:
            continue

        best_by_candidate = {}
        for dist, rec in matches:
            candidate_key = (rec["candidate_type"], rec["candidate_id"])
            current = best_by_candidate.get(candidate_key)
            if current is None or _match_sort_key((dist, rec)) < _match_sort_key(current):
                best_by_candidate[candidate_key] = (dist, rec)

        for candidate_key, (dist, rec) in best_by_candidate.items():
            stats = candidate_stats[candidate_key]
            stats["candidate_type"] = rec["candidate_type"]
            stats["candidate_id"] = rec["candidate_id"]
            stats["matched_reads"] += reads
            stats["exact_reads" if dist == 0 else "hamming1_reads"] += reads
            stats["matched_observed_barcodes"].add(obs_bc)
            stats["matched_expected_barcodes"].add(rec["barcode"])
            match_rows.append({
                "candidate_type": rec["candidate_type"],
                "candidate_id": rec["candidate_id"],
                "wellID": rec["wellID"],
                "barcode_type": rec["barcode_type"],
                "observed_barcode": obs_bc,
                "expected_barcode": rec["barcode"],
                "hamming": dist,
                "reads": reads,
            })

    summaries = []
    for key, stats in candidate_stats.items():
        summaries.append({
            "candidate_type": stats["candidate_type"],
            "candidate_id": stats["candidate_id"],
            "matched_reads": stats["matched_reads"],
            "exact_reads": stats["exact_reads"],
            "hamming1_reads": stats["hamming1_reads"],
            "matched_observed_barcodes": len(stats["matched_observed_barcodes"]),
            "matched_expected_barcodes": len(stats["matched_expected_barcodes"]),
        })
    summaries.sort(key=lambda x: (-x["matched_reads"], -x["matched_expected_barcodes"], _candidate_type_rank(x["candidate_type"]), x["candidate_id"]))

    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    _write_discovery_report(out_file, summaries, match_rows)

    if summaries:
        best = summaries[0]
        top_ties = [
            row for row in summaries
            if row["matched_reads"] == best["matched_reads"]
            and row["matched_expected_barcodes"] == best["matched_expected_barcodes"]
        ]
        tied_types = {row["candidate_type"] for row in top_ties}
        if len(tied_types) > 1:
            labels = ", ".join(
                f"{row['candidate_type']}:{row['candidate_id']}"
                for row in top_ties
            )
            raise ValueError(
                "Barcode discovery could not distinguish manual vs auto candidates. "
                f"Tied candidates: {labels}. Use --manual or --plate to disambiguate. "
                f"See discovery report: {out_file}"
            )

    selected = _select_candidates(summaries, min_unique_barcodes, min_fraction)
    if not selected:
        raise ValueError(
            "Barcode discovery found no confident candidate. "
            f"See discovery report: {out_file}"
        )

    selected_keys = {(row["candidate_type"], row["candidate_id"]) for row in selected}
    selected_records = [
        rec for rec in records
        if (rec["candidate_type"], rec["candidate_id"]) in selected_keys
    ]
    return selected, selected_records


def _manual_records(script_dir, sample_ids=None):
    path = os.path.join(script_dir, "yaml", "manual_barcode_list.yaml")
    with open(path, "r", encoding="utf-8") as handle:
        barcode_set = yaml.safe_load(handle)

    sample_ids = _normalize_ids(sample_ids) if sample_ids is not None else sorted(barcode_set.keys(), key=_natural_key)
    records = []
    for sample_id in sample_ids:
        if sample_id not in barcode_set:
            raise ValueError(f"Manual barcode ID {sample_id} not found in manual_barcode_list.yaml")
        bc_t_i5, bc_n_i5, bc_n_i7 = barcode_set[sample_id]
        for umi_i5, int_i5 in zip(bc_t_i5, bc_n_i5):
            for i7 in bc_n_i7:
                records.append(_record("manual", sample_id, f"MANUAL{sample_id}", umi_i5 + i7, "umi"))
                records.append(_record("manual", sample_id, f"MANUAL{sample_id}", int_i5 + i7, "internal"))
    return records


def _auto_records(script_dir, sample_ids=None):
    path = os.path.join(script_dir, "yaml", "auto_barcode_list.yaml")
    with open(path, "r", encoding="utf-8") as handle:
        barcode_set = yaml.safe_load(handle)

    if sample_ids is None:
        plate_keys = sorted(barcode_set.keys(), key=_natural_key)
    else:
        plate_keys = [f"plate{x}" for x in _normalize_ids(sample_ids)]

    records = []
    for plate_key in plate_keys:
        if plate_key not in barcode_set:
            raise ValueError(f"Plate ID {plate_key.replace('plate', '')} not found in auto_barcode_list.yaml")
        plate_id = plate_key.replace("plate", "")
        for well in sorted(barcode_set[plate_key].keys(), key=_natural_key):
            values = barcode_set[plate_key][well]
            umi_part = values[0] if isinstance(values[0], list) else [values[0]]
            int_part = values[1] if isinstance(values[1], list) else [values[1]]
            well_id = f"P{plate_id}{well}"
            for barcode in umi_part:
                records.append(_record("auto", plate_id, well_id, barcode, "umi"))
            for barcode in int_part:
                records.append(_record("auto", plate_id, well_id, barcode, "internal"))
    return records


def _record(candidate_type, candidate_id, well_id, barcode, barcode_type):
    return {
        "candidate_type": candidate_type,
        "candidate_id": str(candidate_id),
        "wellID": well_id,
        "barcode": str(barcode).strip().upper(),
        "barcode_type": barcode_type,
    }


def _read_bcstats(path):
    counts = {}
    with open(path) as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                count = int(parts[1])
            except ValueError:
                continue
            barcode = parts[0].strip().upper()
            counts[barcode] = counts.get(barcode, 0) + count
    return counts


def _match_sort_key(item):
    dist, rec = item
    return (
        dist,
        rec["wellID"],
        rec["barcode_type"],
        rec["barcode"],
    )


def _candidate_type_rank(candidate_type):
    return {"manual": 0, "auto": 1}.get(candidate_type, 2)


def _write_discovery_report(out_file, summaries, match_rows):
    with open(out_file, "w") as handle:
        print("[summary]", file=handle)
        print("candidate_type\tcandidate_id\tmatched_reads\texact_reads\thamming1_reads\tmatched_observed_barcodes\tmatched_expected_barcodes", file=handle)
        for row in summaries:
            print(
                "\t".join(str(row[col]) for col in (
                    "candidate_type",
                    "candidate_id",
                    "matched_reads",
                    "exact_reads",
                    "hamming1_reads",
                    "matched_observed_barcodes",
                    "matched_expected_barcodes",
                )),
                file=handle,
            )
        print("\n[matches]", file=handle)
        print("candidate_type\tcandidate_id\twellID\tbarcode_type\tobserved_barcode\texpected_barcode\thamming\treads", file=handle)
        for row in sorted(match_rows, key=lambda x: (-x["reads"], x["candidate_type"], x["candidate_id"], x["wellID"])):
            print(
                "\t".join(str(row[col]) for col in (
                    "candidate_type",
                    "candidate_id",
                    "wellID",
                    "barcode_type",
                    "observed_barcode",
                    "expected_barcode",
                    "hamming",
                    "reads",
                )),
                file=handle,
            )


def _select_candidates(summaries, min_unique_barcodes, min_fraction):
    confident = [row for row in summaries if row["matched_expected_barcodes"] >= min_unique_barcodes]
    if not confident:
        return []

    top_type = confident[0]["candidate_type"]
    same_type = [row for row in confident if row["candidate_type"] == top_type]
    best_reads = same_type[0]["matched_reads"]
    threshold = max(1, int(best_reads * min_fraction))
    return [row for row in same_type if row["matched_reads"] >= threshold]


def _normalize_ids(value):
    if isinstance(value, str):
        items = value.split(",")
    else:
        items = value
    return [str(item).strip() for item in items if str(item).strip()]


def _unique(values):
    return list(collections.OrderedDict((v, None) for v in values).keys())


def _natural_key(value):
    text = str(value)
    prefix = "".join(ch for ch in text if not ch.isdigit())
    suffix = "".join(ch for ch in text if ch.isdigit())
    return (prefix, int(suffix) if suffix else -1, text)
