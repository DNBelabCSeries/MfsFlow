"""
Barcode correction utilities for standalone script execution.

This module provides functions for correcting barcodes using barcode maps
and ID maps, enabling accurate cell identification from sequencing data.
"""

from collections import namedtuple
import sys


BarcodeCorrection = namedtuple(
    "BarcodeCorrection",
    ["raw_bc", "corrected_bc", "well_id", "is_internal"],
)


def load_bc_map(binmap_file):
    bc_map = {}
    try:
        with open(binmap_file, "r") as f:
            for line in f:
                if "falseBC" in line:
                    continue
                parts = line.replace(",", "\t").split("\t")
                if len(parts) >= 3:
                    raw = parts[0].strip().upper()
                    fixed = parts[2].strip().upper()
                    if raw and fixed:
                        bc_map[raw] = fixed
    except Exception as e:
        sys.stderr.write(f"Error loading BC map: {e}\n")
    return bc_map


def load_id_map(id_map_file, strict=False):
    id_map = {}
    internal_bcs = set()
    try:
        with open(id_map_file, "r") as f:
            for line in f:
                if line.startswith("wellID"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue

                well_id = parts[0]
                umi_seqs = parts[1].split(",")
                int_seqs = parts[2].split(",")

                for umi_bc in umi_seqs:
                    umi_bc = umi_bc.strip().upper()
                    if umi_bc:
                        id_map[umi_bc] = well_id

                for internal_bc in int_seqs:
                    internal_bc = internal_bc.strip().upper()
                    if internal_bc:
                        id_map[internal_bc] = well_id
                        internal_bcs.add(internal_bc)
    except Exception as e:
        sys.stderr.write(f"Error loading ID map: {e}\n")
        if strict:
            raise
    return id_map, internal_bcs


def correct_read_barcode(read, bc_map, id_map, internal_bcs):
    try:
        value = read.get_tag("CR")
    except KeyError:
        value = None
    raw_bc = value.upper() if isinstance(value, str) else None

    if not raw_bc:
        return None

    corrected_bc = bc_map.get(raw_bc, raw_bc)
    is_internal = corrected_bc in internal_bcs
    _adjust_read_sequence(read, is_internal)

    read.set_tag("CR", raw_bc)
    well_id = id_map.get(corrected_bc) if id_map else None

    if well_id:
        read.set_tag("CC", corrected_bc)
        read.set_tag("CB", well_id)
    else:
        for tag in ("CC", "CB"):
            read.set_tag(tag, None)

    return BarcodeCorrection(
        raw_bc=raw_bc,
        corrected_bc=corrected_bc,
        well_id=well_id,
        is_internal=is_internal,
    )


def _adjust_read_sequence(read, is_internal):
    flag = read.flag
    seq = read.query_sequence
    qual = read.query_qualities

    if flag == 77 and seq:
        if not is_internal:
            if len(seq) > 3:
                read.query_sequence = seq[3:]
                read.query_qualities = qual[3:] if qual is not None and len(qual) == len(seq) else None
            return

        try:
            umi_seq = read.get_tag("UR")
        except KeyError:
            umi_seq = None
        try:
            umi_qual = read.get_tag("UY")
        except KeyError:
            umi_qual = None
        if umi_seq:
            read.query_sequence = umi_seq + seq
            if qual is not None and umi_qual is not None and len(umi_qual) == len(umi_seq):
                read.query_qualities = [ord(c) - 33 for c in umi_qual] + list(qual)
            else:
                read.query_qualities = None
            read.set_tag("UR", None)
            read.set_tag("UY", None)

    elif flag == 141 and is_internal:
        read.set_tag("UR", None)
        read.set_tag("UY", None)
