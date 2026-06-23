#!/usr/bin/env python3
"""
UB tag correction for BAM files in standalone script execution.

This module corrects UMI and barcode tags in BAM files, performing
molecule counting and UMI deduplication for single-cell RNA sequencing
data processing.
"""

import os
import pysam
import argparse
import multiprocessing as mp

_worker_mols = None

def _init_worker(mols_dict):
    global _worker_mols
    _worker_mols = mols_dict

def collect_bam_chunks(inpath, chrs, outpath):
    allpaths = [inpath+".tmp."+c+".bam" for c in chrs[:-1]]
    allpaths.append(inpath+".tmp."+"unmapped"+".bam")
    cat_args = ['-o', outpath]+allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]

def load_bcs(bcpath):
    with open(bcpath) as f:
        x = f.readline()
        y = f.readlines()
        bc = []
        for l in y:
            l = l.split(',')
            bc.append(l[0])
    return(bc)

def load_dict(stub, bcs):
    molecules_dict = {}
    for i in bcs:
        fp = stub+i+".txt"
        if os.path.exists(fp):
            molecules_dict[i] = {}
            with open(fp) as f:
              x = f.readline()
              y = f.readlines()
              for l in y:
                l = l.strip().split('\t')
                if l[3] not in molecules_dict[i]:
                  molecules_dict[i][l[3]] = {}
                if l[0] not in molecules_dict[i][l[3]]:
                  molecules_dict[i][l[3]][l[0]] = {}
                molecules_dict[i][l[3]][l[0]] = l[1]
    return(molecules_dict)

def return_UB(moldict, BC, GE, UX):
    try:
        UB = moldict[BC][GE][UX]
    except KeyError:
        UB = UX
    return(UB)

def correct_tags(args):
    inpath, threads, chr_name = args
    if chr_name == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr_name
    outpath = inpath+".tmp."+chrlabel+".bam"
    with pysam.AlignmentFile(inpath, 'rb', threads=threads) as inp, \
         pysam.AlignmentFile(outpath, 'wb', template=inp, threads=threads) as out:
        for read in inp.fetch(chr_name):
            umi = read.get_tag('UR')
            cell = read.get_tag('CB')
            if read.has_tag('GX'):
                gene = read.get_tag('GX')
            else:
                gene = 'NA'
            umi_new = return_UB(moldict=_worker_mols, BC=cell, GE=gene, UX=umi)
            read.set_tag(tag='UB', value=umi_new, value_type='Z')
            out.write(read)
    return outpath

def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to input BAM file')
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='Path to output bam file')
    parser.add_argument('--p', type=int, default=10,
                        help='Number of processes for bams')
    parser.add_argument('--bcs', type=str, metavar='FILENAME',
                        help='Path to kept barcodes')
    parser.add_argument('--stub', type=str, metavar='FILENAME',
                        help='Molecule table path stub')

    args = parser.parse_args()

    bcs = load_bcs(args.bcs)
    print("Loading molecule correction dictionary...")
    mols = load_dict(args.stub, bcs)
    print("Correcting UB tags...")

    chrs = pysam.idxstats(args.bam).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]

    if args.p > 8:
        pysam_workers = 2
        n_jobs = int(args.p/2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    tasks = [(args.bam, pysam_workers, chr_name) for chr_name in chrs]
    with mp.Pool(n_jobs, initializer=_init_worker, initargs=(mols,)) as pool:
        results = pool.map(correct_tags, tasks)
    _ = results

    collect_bam_chunks(inpath=args.bam, chrs=chrs, outpath=args.out)

if __name__ == "__main__":
    main()
