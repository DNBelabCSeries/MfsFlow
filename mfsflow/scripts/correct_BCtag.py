"""
BC tag correction for BAM files in standalone script execution.

This module corrects barcode tags in BAM files using barcode maps and ID maps,
enabling accurate cell identification and UMI deduplication in single-cell
RNA sequencing data processing.
"""

import sys
import os
try:
    import pysam
except ImportError:
    sys.stderr.write("Error: pysam module is required for this script. Please install it (pip install pysam).\n")
    sys.exit(1)

try:
    from mfsflow.scripts.barcode_corrector import correct_read_barcode, load_bc_map, load_id_map
except ImportError:
    from barcode_corrector import correct_read_barcode, load_bc_map, load_id_map

def main():
    if len(sys.argv) < 6:
        # Args: 1:in_bam, 2:out_bam_umi, 3:out_bam_internal, 4:binmap1, 5+:id_map
        print("Usage: python3 correct_BCtag.py <inbam> <outbam_umi> <outbam_internal> <BCbinmap> [ID_map_file]")
        sys.exit(1)

    in_bam = sys.argv[1]
    out_bam_umi = sys.argv[2]
    out_bam_internal = sys.argv[3]
    binmap = sys.argv[4]  
    id_map_file = sys.argv[5]

    print(f"Loading maps... (BCbinmap: {bool(binmap)}, ID_Map: {bool(id_map_file)})")
    bc_map = load_bc_map(binmap)
    id_map = {}
    internal_bcs = set()
    if id_map_file:
        # strict: a missing/corrupt ID map must abort, otherwise every CB tag
        # would be silently wiped and 100% of reads would end up as Unused BC.
        id_map, internal_bcs = load_id_map(id_map_file, strict=True)
        if not id_map:
            sys.stderr.write(
                f"Error: ID map file loaded but contained no mappings: {id_map_file}\n"
            )
            sys.exit(1)
        print(f"Loaded {len(id_map)} ID mappings and {len(internal_bcs)} internal barcodes.")

    # Open BAM files using context managers for proper resource cleanup
    with pysam.AlignmentFile(in_bam, "rb", check_sq=False) as infile:
        # Prepare header
        header = infile.header.to_dict()
        pg_entry = {
            'ID': 'MfsFlow-correct_BCtag',
            'PN': 'MfsFlow-correct_BCtag',
            'VN': '3.0-pysam-mgi-custom',
            'CL': 'python3 ' + ' '.join(sys.argv)
        }
        if 'PG' in header:
            header['PG'].append(pg_entry)
        else:
            header['PG'] = [pg_entry]

        with pysam.AlignmentFile(out_bam_umi, "wb", header=header) as outfile_umi, \
             pysam.AlignmentFile(out_bam_internal, "wb", header=header) as outfile_int:

            processed = 0

            for read in infile:
                processed += 1
                if processed % 1000000 == 0:
                    print(f"Processed {processed} reads...", flush=True)

                try:
                    correction = correct_read_barcode(read, bc_map, id_map, internal_bcs)
                    if correction is None:
                        outfile_umi.write(read)
                        continue

                    if correction.is_internal:
                        outfile_int.write(read)
                    else:
                        outfile_umi.write(read)

                except Exception as e:
                    # sys.stderr.write(f"Warning: error processing read {read.query_name}: {e}\n")
                    # Write error reads to UMI as safe fallback
                    outfile_umi.write(read)

if __name__ == "__main__":
    main()
