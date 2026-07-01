"""
FASTQ filtering and processing for standalone script execution.

This module provides FASTQ file filtering, quality control, and processing
utilities for single-cell RNA sequencing data, including parallel processing
and quality-based filtering.
"""

import sys
import yaml
import subprocess
import re
import os

import pysam

try:
    from mfsflow.scripts.path_layout import tmp_merge_dir, load_config
except ImportError:
    from path_layout import tmp_merge_dir, load_config

Q30_ASCII = 63
Q30_TABLE = bytes(1 if i >= Q30_ASCII else 0 for i in range(256))

# Alias for backward compatibility
get_config = load_config

def parse_definition(definition):
    if isinstance(definition, list):
        parts = definition
    else:
        parts = str(definition).split(';')
    def_dict = {}
    for part in parts:
        part = str(part).strip()
        match = re.match(r'(\w+)\((.*)\)', part)
        if match:
            key, val = match.groups()
            ranges = []
            for r in val.split(','):
                start, end = map(int, r.split('-'))
                ranges.append((start - 1, end)) # 0-indexed
            def_dict[key] = ranges
    return def_dict

def extract_seq(seq, qual, definition, ss3_no_pattern=False):
    bc_seq = b""
    bc_qual = b""
    umi_seq = b""
    umi_qual = b""
    cdna_seq = b""
    cdna_qual = b""
    
    # Handle BC
    if 'BC' in definition:
        for start, end in definition['BC']:
            bc_seq += seq[start:end]
            bc_qual += qual[start:end]
            
    # Handle UMI
    if 'UMI' in definition:
        if ss3_no_pattern:
            umi_seq = b""
            umi_qual = b""
        else:
            for start, end in definition['UMI']:
                umi_seq += seq[start:end]
                umi_qual += qual[start:end]
                
    # Handle cDNA
    if 'cDNA' in definition:
        # Assuming only one cDNA range as per original perl script logic
        start, end = definition['cDNA'][0]
        if ss3_no_pattern:
            start = 0 # If smart-seq3 pattern not found, take full read
        cdna_seq = seq[start:end]
        cdna_qual = qual[start:end]
        
    return bc_seq, bc_qual, umi_seq, umi_qual, cdna_seq, cdna_qual

def hamming_distance(s1, s2, limit=None):
    if len(s1) != len(s2):
        return len(s1)
    dist = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            dist += 1
            if limit is not None and dist > limit:
                return dist
    return dist

def fastq_iter(handle):
    while True:
        header = handle.readline()
        if not header: break
        seq = handle.readline().rstrip(b'\n\r')
        _plus = handle.readline()
        qual = handle.readline().rstrip(b'\n\r')
        
        if len(seq) != len(qual):
            sys.stderr.write(f"Warning: SEQ/QUAL length mismatch in {header.decode().strip()}: {len(seq)} vs {len(qual)}\n")
            # Pad or truncate qual to match seq length to prevent crashes/offsets
            if len(qual) < len(seq):
                qual += b'#' * (len(seq) - len(qual))
            else:
                qual = qual[:len(seq)]
                
        yield header.rstrip(b'\n\r'), seq, qual

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 fqfilter.py <yaml> <pigz> <tmp_prefix> [--limit N] [--pigz-threads N]")
        sys.exit(1)

    args = sys.argv[1:]
    read_limit = 0
    pigz_threads = 1
    if '--limit' in args:
        try:
            idx = args.index('--limit')
            read_limit = int(args[idx + 1])
            args = args[:idx] + args[idx + 2:]
        except Exception:
            read_limit = 0
    if '--pigz-threads' in args:
        try:
            idx = args.index('--pigz-threads')
            pigz_threads = max(1, int(args[idx + 1]))
            args = args[:idx] + args[idx + 2:]
        except Exception:
            pigz_threads = 1

    yaml_file = args[0]
    pigz = args[1]
    tmp_prefix = args[2]

    config = get_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    
    def sorted_sequence_file_entries(sequence_files):
        if isinstance(sequence_files, dict):
            def key_fn(k):
                m = re.search(r'(\d+)', str(k))
                return int(m.group(1)) if m else str(k)

            for k in sorted(sequence_files.keys(), key=key_fn):
                yield sequence_files[k]
            return
        if isinstance(sequence_files, list):
            for entry in sequence_files:
                yield entry

    seq_entries = list(sorted_sequence_file_entries(config.get('sequence_files', {})))
    filenames = [e['name'] for e in seq_entries if isinstance(e, dict) and e.get('name')]
    base_definitions = [parse_definition(e.get('base_definition', [])) for e in seq_entries if isinstance(e, dict)]
    patterns = [e.get('find_pattern', 'character(0)') for e in seq_entries if isinstance(e, dict)]
    pattern_bytes = []
    for p in patterns:
        if p is None or p == 'character(0)':
            pattern_bytes.append(None)
            continue
        p = str(p)
        mm = None
        if ';' in p:
            p, mm = p.split(';', 1)
            mm = int(mm)
        else:
            mm = 1
        pattern_bytes.append((p.encode('ascii'), mm))

    if isinstance(config.get('filter_cutoffs', {}).get('BC_filter'), dict):
        bc_filter = [
            int(config['filter_cutoffs']['BC_filter']['num_bases']),
            int(config['filter_cutoffs']['BC_filter']['phred']),
        ]
    else:
        bc_filter = list(map(int, str(config['filter_cutoffs']['BC_filter']).split()))

    if isinstance(config.get('filter_cutoffs', {}).get('UMI_filter'), dict):
        umi_filter = [
            int(config['filter_cutoffs']['UMI_filter']['num_bases']),
            int(config['filter_cutoffs']['UMI_filter']['phred']),
        ]
    else:
        umi_filter = list(map(int, str(config['filter_cutoffs']['UMI_filter']).split()))

    def _lowq_table(phred_threshold):
        limit = phred_threshold + 33
        return bytes(1 if i < limit else 0 for i in range(256))

    bc_lowq_table = _lowq_table(bc_filter[1])
    umi_lowq_table = _lowq_table(umi_filter[1])

    merge_dir = tmp_merge_dir(out_dir)
    out_bam = os.path.join(merge_dir, f"{project}{tmp_prefix}.raw.tagged.bam")
    out_bc_stats = os.path.join(merge_dir, f"{project}{tmp_prefix}.BCstats.txt")
    out_q30_stats = os.path.join(merge_dir, f"{project}{tmp_prefix}.Q30stats.txt")

    def split_names(name):
        return [x.strip() for x in str(name).split(',') if x.strip()]

    def make_read_groups():
        groups = []
        fastq_groups = config.get('fastq_groups') or []
        if fastq_groups:
            for row in fastq_groups:
                groups.append({
                    'files': [row['read1'], row['read2']],
                    'fixed_bc': row['barcode'].encode('ascii'),
                })
            return groups

        if len(filenames) < 2:
            groups.append({'files': filenames, 'fixed_bc': None})
            return groups

        r1_files = split_names(filenames[0])
        r2_files = split_names(filenames[1])
        if len(r1_files) != len(r2_files):
            raise ValueError(f"R1/R2 FASTQ count mismatch: {len(r1_files)} vs {len(r2_files)}")
        for r1, r2 in zip(r1_files, r2_files):
            groups.append({'files': [r1, r2], 'fixed_bc': None})
        return groups

    def open_chunk_handle(f, pigz_procs):
        base_name = os.path.basename(f)
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]
        if base_name.endswith('.fastq'):
            base_name = base_name[:-6]
        elif base_name.endswith('.fq'):
            base_name = base_name[:-3]
            
        # Try finding the file with .gz extension first (default behavior)
        chunk_path_gz = os.path.join(merge_dir, f"{base_name}{tmp_prefix}.gz")
        chunk_path_plain = os.path.join(merge_dir, f"{base_name}{tmp_prefix}")
        
        # Determine which file to use
        if os.path.exists(chunk_path_gz):
            proc = subprocess.Popen([pigz, '-p', str(pigz_threads), '-dc', chunk_path_gz], stdout=subprocess.PIPE, text=False, bufsize=1024*1024)
            pigz_procs.append(proc)
            return proc.stdout
        if os.path.exists(chunk_path_plain):
            return open(chunk_path_plain, 'rb')
        return None
    
    bc_stats = {}
    q30_stats = {}

    def add_q30(label, qual):
        if not qual or qual == b"*":
            return
        total_bases, q30_bases = q30_stats.get(label, (0, 0))
        q30_stats[label] = (
            total_bases + len(qual),
            q30_bases + sum(qual.translate(Q30_TABLE)),
        )
    
    bam_header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{
            "ID": "MfsFlow-fqfilter",
            "PN": "MfsFlow-fqfilter",
            "VN": "3.0",
            "CL": "python3 fqfilter.py " + " ".join(sys.argv[1:]),
        }],
    }
    bam_out = pysam.AlignmentFile(out_bam, "wb", header=bam_header)

    total = 0
    filtered = 0
    
    processing_error = None

    def check_qual(q_str, threshold_count, table):
        return sum(q_str.translate(table)) < threshold_count

    def process_records(records, fixed_bc):
        nonlocal total, filtered
        if read_limit > 0 and total >= read_limit:
            return False

        total += 1

        final_bc = b""
        final_bc_q = b""
        final_umi = b""
        final_umi_q = b""
        final_cdna1 = b""
        final_cdna1_q = b""
        final_cdna2 = b""
        final_cdna2_q = b""

        go_ahead = True
        layout = "SE"

        for i, (_header, seq, qual) in enumerate(records):
            read_label = f"R{i + 1}"
            add_q30(read_label, qual)

            ss3_status = "yespattern"
            pat = pattern_bytes[i] if i < len(pattern_bytes) else None
            if pat is not None:
                pat_seq, mm = pat
                if pat_seq == b"ATTGCGCAATG":
                    if hamming_distance(seq[:len(pat_seq)], pat_seq, limit=mm) <= mm:
                        ss3_status = "yespattern"
                    else:
                        ss3_status = "nopattern"
                else:
                    if not seq.startswith(pat_seq):
                        go_ahead = False

            bc, bc_q, umi, umi_q, c1, c1_q = extract_seq(seq, qual, base_definitions[i], ss3_status == "nopattern")

            add_q30(f"{read_label}_BC", bc_q)
            add_q30(f"{read_label}_UMI", umi_q)
            add_q30(f"{read_label}_cDNA", c1_q)
            add_q30("BC", bc_q)
            add_q30("UMI", umi_q)
            add_q30("cDNA", c1_q)

            final_bc += bc
            final_bc_q += bc_q
            final_umi += umi
            final_umi_q += umi_q

            if i == 0:
                final_cdna1, final_cdna1_q = c1, c1_q
            else:
                final_cdna2, final_cdna2_q = c1, c1_q
                layout = "PE"

        if not go_ahead:
            return True

        if fixed_bc:
            final_bc = fixed_bc
            final_bc_q = b"I" * len(fixed_bc)

        if not check_qual(final_bc_q, bc_filter[0], bc_lowq_table):
            return True
        if not check_qual(final_umi_q, umi_filter[0], umi_lowq_table):
            return True

        filtered += 1
        bc_stats[final_bc] = bc_stats.get(final_bc, 0) + 1

        rid = records[0][0].split()[0]
        if rid.startswith(b'@'):
            rid = rid[1:]
        qname = rid.decode('ascii')

        bc_str = final_bc.decode('ascii')
        bc_q_str = final_bc_q.decode('ascii')
        umi_str = final_umi.decode('ascii')
        umi_q_str = final_umi_q.decode('ascii')

        def _make_read(flag, seq, qual):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.flag = flag
            if seq:
                a.query_sequence = seq.decode('ascii')
                a.query_qualities = [b - 33 for b in qual]
            a.set_tag("CR", bc_str, "Z")
            a.set_tag("UR", umi_str, "Z")
            a.set_tag("CY", bc_q_str, "Z")
            a.set_tag("UY", umi_q_str, "Z")
            return a

        if layout == "SE":
            bam_out.write(_make_read(4, final_cdna1, final_cdna1_q))
        else:
            bam_out.write(_make_read(77, final_cdna1, final_cdna1_q))
            bam_out.write(_make_read(141, final_cdna2, final_cdna2_q))
        return True

    try:
        for group in make_read_groups():
            if read_limit > 0 and total >= read_limit:
                break
            pigz_procs = []
            handles = []
            try:
                for f in group['files']:
                    handle = open_chunk_handle(f, pigz_procs)
                    if handle is None:
                        handles = []
                        break
                    handles.append(handle)

                if not handles:
                    for p in pigz_procs:
                        p.wait()
                    continue

                iters = [fastq_iter(h) for h in handles]
                for records in zip(*iters):
                    if not process_records(records, group.get('fixed_bc')):
                        break
            finally:
                for h in handles:
                    try:
                        h.close()
                    except Exception:
                        pass

                for p in pigz_procs:
                    p.wait()
                    if p.returncode not in (0, -13, 141):
                        raise RuntimeError(f"pigz failed (rc={p.returncode}) while reading chunk(s) for prefix {tmp_prefix}")

    except Exception as e:
        processing_error = e
        sys.stderr.write(f"Error: {e}\n")
    finally:
        try:
            bam_out.close()
        except Exception:
            pass

    if processing_error is not None:
        raise processing_error

    # Write BC stats
    with open(out_bc_stats, 'w') as f:
        for bc, count in bc_stats.items():
            f.write(f"{bc.decode('ascii')}\t{count}\n")

    with open(out_q30_stats, 'w') as f:
        f.write("metric\ttotal_bases\tq30_bases\n")
        for metric in sorted(q30_stats):
            total_bases, q30_bases = q30_stats[metric]
            f.write(f"{metric}\t{total_bases}\t{q30_bases}\n")

if __name__ == "__main__":
    main()
