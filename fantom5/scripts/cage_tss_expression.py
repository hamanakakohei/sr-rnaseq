#!/usr/bin/env python3
"""
Compute per-sample CAGE TPM summed over peaks overlapping TSS regions.

Usage:
    python cage_tss_expression.py [options]

Options:
    --gene-list    Gene IDs file (one per line, or TSV with gene ID in last col)
    --gtf          GTF file
    --cage         CAGE peaks TPM file (.osc.txt or .osc.txt.gz)
    --cage-bed     CAGE peaks BED file with hg38 coordinates (col4 = peak name matching --cage)
    --margin       bp margin around TSS (default: 500)
    --out          Output TSV file (default: cage_tss_expression.tsv)
    --bedtools     Path to bedtools binary (default: bedtools)
"""

import argparse
import gzip
import subprocess
import sys
import tempfile
import os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene-list", default="gene_ids.list")
    p.add_argument("--gtf", default="aaa.gtf")
    p.add_argument("--cage", default="hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz")
    p.add_argument("--cage-bed", default="hg38_fair+new_CAGE_peaks_phase1and2.bed.gz")
    p.add_argument("--margin", type=int, default=500)
    p.add_argument("--out", default="cage_tss_expression.tsv")
    p.add_argument("--bedtools", default="bedtools")
    return p.parse_args()


def load_gene_list(path):
    genes = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            # last non-empty column as gene ID
            genes.append(parts[-1].strip())
    return list(dict.fromkeys(genes))  # deduplicated, order preserved


def extract_tss_from_gtf(gtf_path, gene_ids):
    """
    Return dict: gene_id -> list of (chrom, tss_pos, strand)
    TSS = start of transcript (5' end, strand-aware).
    Uses 'transcript' feature lines only.
    """
    gene_set = set(gene_ids)
    tss_dict = defaultdict(set)

    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    with open_fn(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "transcript":
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]

            gene_id = None
            for tok in attrs.split(";"):
                tok = tok.strip()
                if tok.startswith('gene_id'):
                    gene_id = tok.split('"')[1]
                    break
            if gene_id not in gene_set:
                continue

            # GTF is 1-based inclusive; TSS = start for +, end for -
            if strand == "+":
                tss = start - 1  # convert to 0-based
            else:
                tss = end  # 0-based position just after end = 1-based end
            tss_dict[gene_id].add((chrom, tss, strand))

    return tss_dict


def write_tss_bed(tss_dict, margin, path):
    """Write TSS ± margin as BED6 with gene_id in name field."""
    rows = []
    for gene_id, tss_set in tss_dict.items():
        for chrom, tss, strand in tss_set:
            start = max(0, tss - margin)
            end = tss + margin
            rows.append((chrom, start, end, gene_id, ".", strand))

    # Sort for bedtools
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def write_cage_bed_from_bed(cage_bed_path, bed_out):
    """
    Read the hg38 BED file (col0-2: hg38 coords, col3: peak_id, col5: strand).
    Write BED6 for bedtools intersect.
    BED file already uses 0-based half-open coordinates.
    """
    open_fn = gzip.open if cage_bed_path.endswith(".gz") else open
    with open_fn(cage_bed_path, "rt") as f_in, open(bed_out, "w") as f_out:
        for line in f_in:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], parts[1], parts[2]
            peak_id = parts[3]
            strand = parts[5] if len(parts) > 5 else "."
            f_out.write(f"{chrom}\t{start}\t{end}\t{peak_id}\t.\t{strand}\n")


def extract_tpm_from_osc(cage_path, tpm_out):
    """
    Parse CAGE osc file, skip header/stat rows, write peak_id + TPM values.
    Returns list of sample names.
    """
    open_fn = gzip.open if cage_path.endswith(".gz") else open
    sample_names = None
    with open_fn(cage_path, "rt") as f_in, open(tpm_out, "w") as f_out:
        for line in f_in:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            peak_id = parts[0]
            if peak_id == "00Annotation":
                sample_names = parts[1:]
                f_out.write(line)
                continue
            if not (peak_id.startswith("hg19::") or peak_id.startswith("hg38::")):
                continue
            f_out.write(line)
    return sample_names


def run_bedtools_intersect(tss_bed, cage_bed, bedtools):
    """
    Return list of (gene_id, peak_id) for overlapping pairs.
    Uses bedtools intersect -s (strand-aware).
    """
    cmd = [bedtools, "intersect", "-a", tss_bed, "-b", cage_bed,
           "-s", "-wa", "-wb"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"bedtools error:\n{result.stderr}")

    pairs = []
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        gene_id = cols[3]
        peak_id = cols[9]  # name field of cage_bed (col 4 = index 3, but -wb gives all cage cols after tss cols)
        pairs.append((gene_id, peak_id))
    return pairs


def load_tpm_for_peaks(tpm_file, peak_ids):
    """Return dict peak_id -> list[float] and sample names."""
    wanted = set(peak_ids)
    tpm_data = {}
    sample_names = None
    with open(tpm_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "00Annotation":
                sample_names = parts[1:]
                continue
            if parts[0] in wanted:
                tpm_data[parts[0]] = [float(v) for v in parts[1:]]
    return sample_names, tpm_data


def main():
    args = parse_args()

    print("Loading gene list...", flush=True)
    gene_ids = load_gene_list(args.gene_list)
    print(f"  {len(gene_ids)} genes", flush=True)

    print("Extracting TSS from GTF...", flush=True)
    tss_dict = extract_tss_from_gtf(args.gtf, gene_ids)
    missing = [g for g in gene_ids if g not in tss_dict]
    if missing:
        print(f"  WARNING: {len(missing)} genes not found in GTF: {missing[:5]}{'...' if len(missing)>5 else ''}", flush=True)
    print(f"  {len(tss_dict)} genes with TSS", flush=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tss_bed = os.path.join(tmpdir, "tss.bed")
        cage_bed = os.path.join(tmpdir, "cage.bed")
        cage_tpm = os.path.join(tmpdir, "cage_tpm.tsv")

        print(f"Writing TSS BED (margin ±{args.margin} bp)...", flush=True)
        write_tss_bed(tss_dict, args.margin, tss_bed)

        print("Converting CAGE BED to hg38 coordinates...", flush=True)
        write_cage_bed_from_bed(args.cage_bed, cage_bed)
        print(f"  {sum(1 for _ in open(cage_bed))} peaks", flush=True)

        print("Parsing CAGE TPM file...", flush=True)
        sample_names = extract_tpm_from_osc(args.cage, cage_tpm)
        print(f"  {len(sample_names)} samples", flush=True)

        print("Running bedtools intersect...", flush=True)
        pairs = run_bedtools_intersect(tss_bed, cage_bed, args.bedtools)
        print(f"  {len(pairs)} TSS-peak overlaps", flush=True)

        if not pairs:
            sys.exit("No overlaps found. Check coordinates or increase --margin.")

        # gene_id -> set of peak_ids
        gene_peaks = defaultdict(set)
        for gene_id, peak_id in pairs:
            gene_peaks[gene_id].add(peak_id)

        all_peak_ids = set(pid for pids in gene_peaks.values() for pid in pids)

        print("Loading TPM values...", flush=True)
        sample_names, tpm_data = load_tpm_for_peaks(cage_tpm, all_peak_ids)

        print(f"Writing output to {args.out}...", flush=True)
        n_samples = len(sample_names)
        with open(args.out, "w") as f:
            f.write("gene_id\toverlapping_peaks\t" + "\t".join(sample_names) + "\n")
            for gene_id in gene_ids:
                if gene_id not in gene_peaks:
                    f.write(gene_id + "\tNA\t" + "\t".join(["0"] * n_samples) + "\n")
                    continue
                peaks = sorted(gene_peaks[gene_id])
                # Sum TPM across overlapping peaks
                totals = [0.0] * n_samples
                for pid in peaks:
                    if pid in tpm_data:
                        for i, v in enumerate(tpm_data[pid]):
                            totals[i] += v
                f.write(gene_id + "\t" + ",".join(peaks) + "\t"
                        + "\t".join(f"{v:.6g}" for v in totals) + "\n")

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
