"""
Stream correction for BAM files in standalone script execution.

This module provides streaming barcode and UMI correction for BAM files,
enabling real-time correction of sequencing reads during processing for
single-cell RNA sequencing data.
"""

import sys
import os
try:
    import pysam
except ImportError:
    pysam = None

try:
    from mfsflow.scripts.barcode_corrector import BarcodeCorrection, correct_read_barcode, load_bc_map, load_id_map
except ImportError:
    from barcode_corrector import BarcodeCorrection, correct_read_barcode, load_bc_map, load_id_map


def get_or_apply_correction(read, bc_map, id_map, internal_bcs):
    if read.has_tag("CC") and read.has_tag("CB"):
        corrected_bc = read.get_tag("CC")
        raw_bc = read.get_tag("CR") if read.has_tag("CR") else corrected_bc
        return BarcodeCorrection(
            raw_bc=raw_bc,
            corrected_bc=corrected_bc,
            well_id=read.get_tag("CB"),
            is_internal=corrected_bc in internal_bcs,
        )
    return correct_read_barcode(read, bc_map, id_map, internal_bcs)

def main():
    if pysam is None:
        sys.stderr.write("Error: pysam module is required for this script. Please install it (pip install pysam).\n")
        sys.exit(1)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--binning', required=True, help="Barcode binning file")
    parser.add_argument('--idmap', required=True, help="ID map file")
    parser.add_argument('--type', choices=['umi', 'internal'], required=True, help="Output type filter")
    parser.add_argument('bam_files', nargs='*', default=['-'], help="Input BAM files")
    args = parser.parse_args()

    # Load Maps
    bc_map = load_bc_map(args.binning)
    id_map, internal_bcs = load_id_map(args.idmap, strict=True)
    target_type = args.type

    # Open Output (Standard Output) - initialized on first file
    outfile = None

    try:
        for bam_path in args.bam_files:
            # Handle '-' for stdin
            if bam_path == '-':
                f_obj = sys.stdin.buffer
            else:
                f_obj = bam_path # pysam accepts path string

            try:
                infile = pysam.AlignmentFile(f_obj, "rb", check_sq=False)
            except ValueError:
                continue

            # Initialize output using header from first file
            # Output BAM (binary) to stdout; STAR reads it via --readFilesCommand samtools view.
            # Binary BAM avoids the ~30% SAM text serialization overhead of pysam mode "w".
            if outfile is None:
                try:
                    # Mode "wb" = BAM binary. File "-" = stdout.
                    outfile = pysam.AlignmentFile("-", "wb", template=infile)
                except (BrokenPipeError, IOError):
                    infile.close()
                    return

            try:
                for read in infile:
                    try:
                        correction = get_or_apply_correction(read, bc_map, id_map, internal_bcs)
                        if correction is None:
                            if target_type == 'umi':
                                outfile.write(read)
                            continue

                        # Filter Logic: Only output if matches target type
                        if target_type == 'umi' and correction.is_internal:
                            continue
                        if target_type == 'internal' and not correction.is_internal:
                            continue

                        outfile.write(read)

                    except (BrokenPipeError, IOError) as e:
                        if getattr(e, "errno", None) == 32: # EPIPE
                            break
                        else:
                            raise
                    except Exception:
                        if target_type == 'umi':
                            outfile.write(read)
            
            finally:
                infile.close()

    except (BrokenPipeError, KeyboardInterrupt):
        pass
    finally:
        if outfile:
            try:
                outfile.close()
            except (BrokenPipeError, IOError):
                pass

if __name__ == "__main__":
    main()
