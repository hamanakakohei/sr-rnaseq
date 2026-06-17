#!/usr/bin/env python3
"""
Calculate total BigWig signal over all exon regions for specified transcripts or genes.

Notes:
    - tx mode: sums signal over each transcript's exons
    - gene mode: merges overlapping exons across all transcripts per gene before summing
                 (avoids double-counting shared exon regions)
"""

import argparse
import sys
import pyBigWig
import pyranges as pr


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bw",   required=True, help="Input BigWig file")
    p.add_argument("--gtf",  required=True, help="Input GTF file")
    p.add_argument("--list", required=True, help="File with IDs (one per line)", dest="id_list")
    p.add_argument("--out",  required=True, help="Output TSV file")
    p.add_argument("--mode", choices=["tx", "gene"], default="gene",
                   help="ID type in --list: 'tx' for transcript_id, 'gene' for gene_id")
    return p.parse_args()


def sum_signal(bw, intervals_df):
    total = 0.0
    for _, row in intervals_df.iterrows():
        val = bw.stats(row.Chromosome, int(row.Start), int(row.End), type="sum")[0]
        total += val if val is not None else 0.0
    return total


def main():
    args = parse_args()

    with open(args.id_list) as f:
        ids = [line.strip() for line in f if line.strip()]
    if not ids:
        sys.exit("ERROR: No IDs found in list file.")
    print(f"Loaded {len(ids)} {args.mode} ID(s).", file=sys.stderr)

    gtf = pr.read_gtf(args.gtf)
    exons = gtf[gtf.Feature == "exon"]

    id_col = "transcript_id" if args.mode == "tx" else "gene_id"
    exons = exons[exons.df[id_col].isin(ids)]
    if len(exons) == 0:
        sys.exit(f"ERROR: No exon records found for specified {args.mode} IDs.")

    bw = pyBigWig.open(args.bw)
    results = {}

    for query_id in ids:
        subset = exons[exons.df[id_col] == query_id]
        if len(subset) == 0:
            print(f"WARNING: No exons found for {query_id}", file=sys.stderr)
            results[query_id] = "NA"
            continue

        if args.mode == "gene":
            # 複数転写産物由来の重複エクソンをマージしてから集計
            subset = subset.merge()

        total = sum_signal(bw, subset.df)
        results[query_id] = total
        print(f"  {query_id}: {total}", file=sys.stderr)

    bw.close()

    header = "transcript_id" if args.mode == "tx" else "gene_id"
    with open(args.out, "w") as out:
        out.write(f"{header}\ttotal_signal\n")
        for query_id in ids:
            out.write(f"{query_id}\t{results[query_id]}\n")

    print(f"Done. Results written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
