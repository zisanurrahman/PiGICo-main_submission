#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 23:09:02 2025

@author: Zisanur Rahman

Batch-capable RPKM calculator with .cover-based filtering.

Features:
- Accepts a single BAM(.gz) or a directory of BAM(.gz) files.
- Accepts a single .cover(.gz) or a directory of .cover(.gz) files.
- Matches cover files to BAMs by basename (sample.bam -> sample.cover).
- Filters contigs by coverage, meandepth, covbases from the .cover file.
- Computes RPKM on kept contigs.
- Outputs TSVs with extra columns:
  Contig, Length(bp), Mapped_Reads, RPKM, Coverage, MeanDepth, CovBases

"""


import os
import sys
import argparse
import glob
import gzip
from collections import defaultdict

import pysam


# ------------------------- I/O helpers -------------------------

def open_maybe_gzip(path, mode='rt'):
    """Open plain text or gzipped text file."""
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def get_contig_lengths(fai_path):
    """Parse FASTA .fai index -> dict {contig: length_bp}."""
    contig_lengths = {}
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                contig, length = parts[0], parts[1]
                try:
                    contig_lengths[contig] = int(length)
                except ValueError:
                    pass
    if not contig_lengths:
        sys.exit(f"[ERROR] No contig lengths parsed from {fai_path}")
    return contig_lengths


def count_reads_per_contig(bam_path):
    """Count mapped reads per contig from BAM or BAM.GZ."""
    bam = pysam.AlignmentFile(bam_path, "rb")  # pysam auto-detects .bam or .bam.gz
    counts = defaultdict(int)
    total = 0
    for read in bam.fetch(until_eof=True):
        if not read.is_unmapped:
            counts[read.reference_name] += 1
            total += 1
    bam.close()
    return counts, total


# ------------------------- .cover parsing -------------------------

def load_cover_filter(cover_path, min_coverage, min_meandepth, min_covbases):
    """
    Read .cover(.gz) and return:
      keep_info: dict { contig: (coverage, meandepth, covbases) } for contigs passing thresholds.
    """
    keep_info = {}
    with open_maybe_gzip(cover_path, 'rt') as f:
        header = None
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                header = [h.lstrip('#') for h in line.split('\t')]
                continue

            parts = line.split('\t')
            if header and len(parts) == len(header):
                row = dict(zip(header, parts))
                try:
                    rname     = row.get('rname') or row.get('#rname') or parts[0]
                    coverage  = float(row.get('coverage', parts[5]))
                    meandepth = float(row.get('meandepth', parts[6]))
                    covbases  = float(row.get('covbases', parts[4]))
                except (ValueError, IndexError):
                    continue
            else:
                try:
                    rname     = parts[0]
                    covbases  = float(parts[4])
                    coverage  = float(parts[5])
                    meandepth = float(parts[6])
                except (ValueError, IndexError):
                    continue

            if (coverage >= min_coverage and
                meandepth >= min_meandepth and
                covbases >= min_covbases):
                keep_info[rname] = (coverage, meandepth, covbases)
    return keep_info


# ------------------------- RPKM calculation -------------------------

def compute_rpkm(read_counts, contig_lengths, total_mapped_reads):
    rpkm = {}
    if total_mapped_reads <= 0:
        return rpkm
    for contig, count in read_counts.items():
        L_kb = contig_lengths.get(contig, 0) / 1000.0
        if L_kb > 0:
            rpkm[contig] = (count / (L_kb * total_mapped_reads)) * 1e6
    return rpkm


# ------------------------- Input resolution -------------------------

def resolve_inputs(b_arg, c_arg):
    """Return list of (bam_path, cover_path, sample) tuples."""
    pairs = []
    if os.path.isdir(b_arg):
        bam_list = sorted(glob.glob(os.path.join(b_arg, '*.bam')) +
                          glob.glob(os.path.join(b_arg, '*.bam.gz')))
        if not bam_list:
            sys.exit(f"No BAMs found in directory: {b_arg}")
        for bam in bam_list:
            sample = os.path.splitext(os.path.basename(bam))[0]
            if sample.endswith('.bam'):
                sample = sample[:-4]
            # Find matching cover
            if os.path.isdir(c_arg):
                cov = os.path.join(c_arg, f'{sample}.cover')
                if not os.path.exists(cov):
                    cov_gz = cov + '.gz'
                    cov = cov_gz if os.path.exists(cov_gz) else None
            else:
                cov = c_arg if os.path.isfile(c_arg) else None
            if not cov:
                print(f"[WARN] Cover missing for {sample}; skipping", file=sys.stderr)
                continue
            pairs.append((bam, cov, sample))
    else:
        if not os.path.isfile(b_arg):
            sys.exit(f"BAM not found: {b_arg}")
        if not os.path.isfile(c_arg):
            sys.exit(f"COVER not found: {c_arg}")
        sample = os.path.basename(b_arg)
        for ext in ('.bam.gz', '.bam'):
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        pairs.append((b_arg, c_arg, sample))
    return pairs


# ------------------------- Main -------------------------

def main():
    ap = argparse.ArgumentParser(description="Calculate RPKM with .cover-based filtering (single sample or batch).")
    ap.add_argument("-b", "--bam", required=True, help="BAM(.gz) file OR directory of BAM(.gz) files")
    ap.add_argument("-c", "--cover", required=True, help=".cover(.gz) file OR directory of .cover(.gz) files")
    ap.add_argument("-f", "--fai", required=True, help="Reference FASTA .fai")
    ap.add_argument("-o", "--out", required=True, help="Output file (single) OR output directory (batch)")
    ap.add_argument("--min-coverage", type=float, default=80.0, help="Min coverage %% (default 80)")
    ap.add_argument("--min-meandepth", type=float, default=1.0, help="Min mean depth (default 1)")
    ap.add_argument("--min-covbases", type=float, default=500.0, help="Min covbases (default 500)")
    ap.add_argument("--denominator", choices=["filtered","all"], default="filtered",
                    help="Normalize by reads on filtered contigs (filtered) or all contigs (all)")
    args = ap.parse_args()

    contig_lengths = get_contig_lengths(args.fai)
    batch = os.path.isdir(args.bam)

    pairs = resolve_inputs(args.bam, args.cover)

    if batch:
        out_dir = args.out
        os.makedirs(out_dir, exist_ok=True)
    else:
        out_file = args.out

    for bam_path, cover_path, sample in pairs:
        print(f"▶ Processing {sample}")
        read_counts_all, total_mapped_all = count_reads_per_contig(bam_path)
        keep_info = load_cover_filter(cover_path, args.min_coverage, args.min_meandepth, args.min_covbases)
        read_counts_kept = {c: n for c, n in read_counts_all.items() if c in keep_info}

        total = sum(read_counts_kept.values()) if args.denominator == "filtered" else total_mapped_all
        rpkm = compute_rpkm(read_counts_kept, contig_lengths, total)

        if batch:
            out_path = os.path.join(args.out, f"{sample}.rpkm.tsv")
        else:
            out_path = out_file

        with open(out_path, 'w') as out:
            out.write("# thresholds: coverage>={:.2f} meandepth>={:.2f} covbases>={:.0f}; denominator={}\n"
                      .format(args.min_coverage, args.min_meandepth, args.min_covbases, args.denominator))
            out.write("Contig\tLength(bp)\tMapped_Reads\tRPKM\tCoverage\tMeanDepth\tCovBases\n")
            for contig in sorted(rpkm):
                length_bp = contig_lengths.get(contig, 0)
                mapped    = read_counts_kept.get(contig, 0)
                cov, depth, bases = keep_info.get(contig, (0,0,0))
                out.write(f"{contig}\t{length_bp}\t{mapped}\t{rpkm[contig]:.4f}\t{cov:.2f}\t{depth:.2f}\t{bases:.0f}\n")
        print(f"✅ Done: {out_path}")


if __name__ == "__main__":
    main()

