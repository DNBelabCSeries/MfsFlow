"""
Stream correction for BAM files in standalone script execution.

This module provides streaming barcode and UMI correction for BAM files,
enabling real-time correction of sequencing reads during processing for
single-cell RNA sequencing data.
"""

import sys
try:
    import pysam
except ImportError:
    pysam = None

try:
    from mfsflow.scripts.barcode_corrector import BarcodeCorrection, correct_read_barcode, load_bc_map, load_id_map
except ImportError:
    from barcode_corrector import BarcodeCorrection, correct_read_barcode, load_bc_map, load_id_map


def get_or_apply_correction(read, bc_map, id_map, internal_bcs):
    # Fetch each tag once. This function runs for every BAM record in the
    # streaming Mapping path, so has_tag()+get_tag() pairs add measurable C
    # extension overhead on large libraries.
    try:
        corrected_bc = read.get_tag("CC")
        well_id = read.get_tag("CB")
    except KeyError:
        return correct_read_barcode(read, bc_map, id_map, internal_bcs)
    try:
        raw_bc = read.get_tag("CR")
    except KeyError:
        raw_bc = corrected_bc
    if corrected_bc and well_id:
        return BarcodeCorrection(
            raw_bc=raw_bc,
            corrected_bc=corrected_bc,
            well_id=well_id,
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

    input_reads = 0
    output_reads = 0
    filtered_reads = 0
    broken_pipe = False

    try:
        for bam_path in args.bam_files:
            # Handle '-' for stdin
            if bam_path == '-':
                f_obj = sys.stdin.buffer
            else:
                f_obj = bam_path # pysam accepts path string

            try:
                infile = pysam.AlignmentFile(f_obj, "rb", check_sq=False)
            except Exception as exc:
                raise RuntimeError(f"Unable to open input BAM {bam_path}: {exc}") from exc

            # Initialize output using header from first file
            # Output BAM (binary) to stdout; STAR reads it via --readFilesCommand samtools view.
            # Binary BAM avoids the ~30% SAM text serialization overhead of pysam mode "w".
            if outfile is None:
                try:
                    # Mode "wb" = BAM binary. File "-" = stdout.
                    outfile = pysam.AlignmentFile("-", "wb", template=infile)
                except (BrokenPipeError, IOError) as exc:
                    infile.close()
                    if isinstance(exc, BrokenPipeError) or getattr(exc, "errno", None) == 32:
                        return 141
                    raise

            try:
                for read in infile:
                    input_reads += 1
                    try:
                        correction = get_or_apply_correction(read, bc_map, id_map, internal_bcs)
                        if correction is None:
                            if target_type == 'umi':
                                outfile.write(read)
                                output_reads += 1
                            else:
                                filtered_reads += 1
                            continue

                        # Filter Logic: Only output if matches target type
                        if target_type == 'umi' and correction.is_internal:
                            filtered_reads += 1
                            continue
                        if target_type == 'internal' and not correction.is_internal:
                            filtered_reads += 1
                            continue

                        outfile.write(read)
                        output_reads += 1

                    except (BrokenPipeError, IOError) as e:
                        if isinstance(e, BrokenPipeError) or getattr(e, "errno", None) == 32:
                            raise BrokenPipeError("STAR closed the correction stream early") from e
                        raise
                    except Exception as exc:
                        read_name = getattr(read, "query_name", "<unknown>")
                        raise RuntimeError(f"Barcode correction failed for read {read_name}: {exc}") from exc
            
            finally:
                infile.close()

    except BrokenPipeError:
        broken_pipe = True
    except KeyboardInterrupt:
        return 130
    finally:
        if outfile:
            try:
                outfile.close()
            except (BrokenPipeError, IOError) as exc:
                if isinstance(exc, BrokenPipeError) or getattr(exc, "errno", None) == 32:
                    broken_pipe = True
                else:
                    raise

    sys.stderr.write(
        f"stream_corrector summary: input={input_reads}, output={output_reads}, "
        f"filtered={filtered_reads}, type={target_type}\n"
    )
    return 141 if broken_pipe else 0

if __name__ == "__main__":
    sys.exit(main())
