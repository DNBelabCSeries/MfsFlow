"""
Pipeline configuration builder: parses CLI arguments, detects read layout, and resolves reference files.

This module constructs the pipeline configuration from command-line arguments,
detects read layout (SE/PE), validates input files, and resolves reference
file paths for the analysis pipeline.
"""

import copy
import csv
from datetime import datetime
import gzip
import glob
import logging
import os

from mfsflow.config.defaults import DEFAULT_CONFIG

logger = logging.getLogger(__name__)


def build_base_config(args, script_dir):
    """Build base pipeline configuration from command-line arguments.
    
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        script_dir (str): Directory containing pipeline scripts.
        
    Returns:
        tuple: (config_dict, samplesheet_records).
    """
    config = copy.deepcopy(DEFAULT_CONFIG)
    config["project"] = args.sample
    config["num_threads"] = args.threads
    config["toolkit_directory"] = script_dir
    config["which_Stage"] = args.stage
    if args.tmpRoot:
        config["performance_opts"]["tmp_root"] = os.path.abspath(args.tmpRoot)

    root_out = os.path.abspath(args.outdir) if args.outdir else os.path.join(os.getcwd(), args.sample)
    config["out_dir"] = os.path.join(root_out, "XPRESS_PROCESSING")

    if args.manual:
        config["sample"]["sample_type"] = "manual"
        config["sample"]["sample_id"] = args.manual
    elif args.plate:
        config["sample"]["sample_type"] = "auto"
        config["sample"]["sample_id"] = args.plate
    elif args.expectBarcode:
        config["sample"]["sample_type"] = "custom"
        config["sample"]["sample_id"] = "1"
        if not os.path.exists(args.expectBarcode):
            raise FileNotFoundError(f"Custom barcode file not found: {args.expectBarcode}")
        config["barcodes"]["barcode_file"] = os.path.abspath(args.expectBarcode)
    else:
        # Default: discover mode — infer plate/manual barcode set from reads
        config["sample"]["sample_type"] = "discover"
        config["sample"]["sample_id"] = "discover"

    samplesheet_records = load_samplesheet(args.samplesheet, args.fastqs)
    if samplesheet_records:
        fastq_pairs = [(r["read1"], r["read2"]) for r in samplesheet_records]
    else:
        fastq_pairs = discover_fastq_pairs(args.fastqs)

    configure_reads(config, fastq_pairs, bool(samplesheet_records))
    configure_reference(config, args.genomeDir)
    return config, samplesheet_records


def configure_reads(config, fastq_pairs, has_samplesheet):
    """Configure read files and detect read layout.
    
    Args:
        config (dict): Pipeline configuration to update.
        fastq_pairs (list): List of (R1, R2) FASTQ file paths.
        has_samplesheet (bool): Whether samplesheet was provided.
    """
    config["sequence_files"]["file1"]["name"] = ",".join(r1 for r1, _r2 in fastq_pairs)
    config["sequence_files"]["file2"]["name"] = ",".join(r2 for _r1, r2 in fastq_pairs)

    logger.info("Detecting read lengths...")
    len_r1 = get_read_length(fastq_pairs[0][0])
    len_r2 = get_read_length(fastq_pairs[0][1])
    logger.info(f"Detected R1 Length: {len_r1}, R2 Length: {len_r2}")

    for r1_file, r2_file in fastq_pairs[1:]:
        cur_r1 = get_read_length(r1_file)
        cur_r2 = get_read_length(r2_file)
        if cur_r1 != len_r1 or cur_r2 != len_r2:
            raise ValueError(
                f"All FASTQ pairs must have identical read lengths. "
                f"Expected R1/R2 {len_r1}/{len_r2}, got {cur_r1}/{cur_r2} for {r1_file}, {r2_file}"
            )

    umi_len = 10
    bc_len = 20
    if len_r1 <= umi_len:
        raise ValueError(f"R1 length ({len_r1}) must be > UMI length ({umi_len}).")

    if len_r2 == len_r1 + bc_len:
        input_mode = "read_embedded_barcode"
    elif len_r2 == len_r1:
        if not has_samplesheet:
            raise ValueError("R1/R2 lengths are equal. Provide --samplesheet so barcode can be assigned per FASTQ pair.")
        input_mode = "samplesheet_barcode"
    else:
        raise ValueError(
            f"Unsupported read length structure: R1={len_r1}, R2={len_r2}. "
            f"Expected R2=R1+20 or R2=R1 with --samplesheet."
        )

    config["barcode_source"] = input_mode
    config["sequence_files"]["file1"]["base_definition"] = [
        f"cDNA({umi_len + 1}-{len_r1})",
        f"UMI(1-{umi_len})",
    ]
    if input_mode == "read_embedded_barcode":
        r2_cdna_end = len_r2 - bc_len
        r2_bc_start = len_r2 - bc_len + 1
        config["sequence_files"]["file2"]["base_definition"] = [
            f"cDNA(1-{r2_cdna_end})",
            f"BC({r2_bc_start}-{len_r2})",
        ]
    else:
        config["sequence_files"]["file2"]["base_definition"] = [
            f"cDNA(1-{len_r2})",
        ]


def configure_reference(config, genome_dir):
    """Configure reference genome paths.
    
    Args:
        config (dict): Pipeline configuration to update.
        genome_dir (str): Directory containing reference genome files.
    """
    genome_dir = os.path.abspath(genome_dir)
    config["reference"]["STAR_index"] = os.path.join(genome_dir, "star")

    gtf_candidates = [
        os.path.join(genome_dir, "genes", "genes.gtf"),
        os.path.join(genome_dir, "genes", "genes.gtf.gz"),
    ]
    gtf_file = next((path for path in gtf_candidates if os.path.exists(path)), None)
    if not gtf_file:
        raise FileNotFoundError(
            "GTF file not found. Checked: " + ", ".join(gtf_candidates)
        )
    config["reference"]["GTF_file"] = gtf_file


def get_read_length(fastq_path):
    """Get the read length from a FASTQ file.
    
    Args:
        fastq_path (str): Path to FASTQ file (optionally gzipped).
        
    Returns:
        int: Read length in bases.
        
    Raises:
        ValueError: If read length cannot be determined.
    """
    try:
        opener = gzip.open if fastq_path.endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            f.readline()
            seq = f.readline().strip()
            if not seq:
                raise ValueError(f"Empty sequence in {fastq_path}")
            return len(seq)
    except Exception as e:
        raise ValueError(f"Failed to determine read length for {fastq_path}: {e}")


def discover_fastq_pairs(fastqs_dir):
    """Discover R1/R2 FASTQ pairs in a directory.
    
    Args:
        fastqs_dir (str): Directory containing FASTQ files.
        
    Returns:
        list: List of (R1_path, R2_path) tuples.
        
    Raises:
        NotADirectoryError: If directory doesn't exist.
        FileNotFoundError: If no pairs found or missing R2 files.
    """
    fastqs_dir = os.path.abspath(fastqs_dir)
    if not os.path.isdir(fastqs_dir):
        raise NotADirectoryError(f"FASTQ directory not found: {fastqs_dir}")

    candidates = []
    for pattern in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"):
        candidates.extend(glob.glob(os.path.join(fastqs_dir, pattern)))
    candidates = sorted(set(candidates))

    r1_files = [p for p in candidates if _is_r1_fastq(os.path.basename(p))]
    pairs = []
    missing_pairs = []
    for r1 in r1_files:
        r2 = _guess_r2_path(r1, candidates)
        if r2:
            pairs.append((os.path.abspath(r1), os.path.abspath(r2)))
        else:
            missing_pairs.append(os.path.abspath(r1))

    if missing_pairs:
        raise FileNotFoundError(
            "Missing matching R2 FASTQ for: " + ", ".join(missing_pairs)
        )

    if not pairs:
        raise FileNotFoundError(f"No R1/R2 FASTQ pairs found in {fastqs_dir}")
    return pairs


def _is_r1_fastq(name):
    tokens = ("_R1", ".R1", "_1.", ".1.")
    return any(t in name for t in tokens)


def _guess_r2_path(r1, candidates):
    name = os.path.basename(r1)
    replacements = (
        ("_R1", "_R2"),
        (".R1", ".R2"),
        ("_1.", "_2."),
        (".1.", ".2."),
    )
    candidate_set = {os.path.abspath(p) for p in candidates}
    for old, new in replacements:
        if old in name:
            guessed = os.path.abspath(os.path.join(os.path.dirname(r1), name.replace(old, new, 1)))
            if guessed in candidate_set:
                return guessed
    return None


def load_samplesheet(samplesheet, fastqs_dir):
    """Load samplesheet CSV and validate FASTQ paths.
    
    Args:
        samplesheet (str): Path to samplesheet CSV file.
        fastqs_dir (str): Directory containing FASTQ files.
        
    Returns:
        list: List of sample record dictionaries.
        
    Raises:
        ValueError: If samplesheet is invalid or has missing records.
    """
    if not samplesheet:
        return []
    samplesheet = os.path.abspath(samplesheet)
    fastqs_dir = os.path.abspath(fastqs_dir)
    records = []
    seen_pairs = set()
    with open(samplesheet, newline="") as f:
        reader = csv.DictReader(f)
        required = {"read1", "read2", "barcode"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Samplesheet missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            read1 = _resolve_samplesheet_path(row["read1"], fastqs_dir)
            read2 = _resolve_samplesheet_path(row["read2"], fastqs_dir)
            barcode = row["barcode"].strip().upper()
            if not barcode:
                raise ValueError(f"Samplesheet barcode must not be empty: {row}")
            pair_key = (read1, read2)
            if pair_key in seen_pairs:
                raise ValueError(f"Duplicate read1/read2 pair in samplesheet: {read1}, {read2}")
            seen_pairs.add(pair_key)
            records.append({
                "read1": read1,
                "read2": read2,
                "barcode": barcode,
            })

    if not records:
        raise ValueError(f"Samplesheet has no records: {samplesheet}")
    return records


def _resolve_samplesheet_path(path, fastqs_dir):
    path = path.strip()
    resolved = path if os.path.isabs(path) else os.path.join(fastqs_dir, path)
    resolved = os.path.abspath(resolved)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"FASTQ listed in samplesheet not found: {resolved}")
    return resolved


def load_expected_barcodes(expect_id_file):
    """Load expected barcode table from TSV file.
    
    Args:
        expect_id_file (str): Path to expected barcode TSV file.
        
    Returns:
        dict: Dictionary mapping well IDs to barcode lists.
    """
    expected = {}
    with open(expect_id_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "wellID":
                continue
            expected[parts[0]] = {
                "umi": [x.strip().upper() for x in parts[1].split(",") if x.strip()],
                "internal": [x.strip().upper() for x in parts[2].split(",") if x.strip()],
            }
    return expected


def resolve_samplesheet_barcodes(records, expect_id_file):
    """Resolve samplesheet barcodes against expected barcode table.
    
    Args:
        records (list): List of sample record dictionaries.
        expect_id_file (str): Path to expected barcode TSV file.
        
    Returns:
        list: Updated records with well ID and barcode type information.
        
    Raises:
        ValueError: If barcode not found or duplicated.
    """
    expected = load_expected_barcodes(expect_id_file)
    barcode_lookup = {}
    for well_id, values in expected.items():
        for barcode_type in ("umi", "internal"):
            for barcode in values[barcode_type]:
                if barcode in barcode_lookup:
                    prev = barcode_lookup[barcode]
                    raise ValueError(
                        f"Barcode {barcode} is duplicated in expected barcode table: "
                        f"{prev['wellID']} and {well_id}"
                    )
                barcode_lookup[barcode] = {
                    "wellID": well_id,
                    "barcode_type": barcode_type,
                }

    resolved = []
    for row in records:
        barcode = row["barcode"]
        if barcode not in barcode_lookup:
            raise ValueError(f"Samplesheet barcode not found in expected barcode table: {barcode}")

        new_row = dict(row)
        new_row.update(barcode_lookup[barcode])
        resolved.append(new_row)
    return resolved


def resolve_samplesheet_fastq_groups(records, expect_id_file, sample_type):
    """Build FASTQ groups for samplesheet mode.

    In discover mode the barcode can legitimately appear in both manual and
    auto candidate tables. If discovery has not already narrowed the expected
    table, keep the user-provided barcode unresolved.
    """
    if str(sample_type or "").strip().lower() == "discover":
        return [dict(row) for row in records]
    return resolve_samplesheet_barcodes(records, expect_id_file)
