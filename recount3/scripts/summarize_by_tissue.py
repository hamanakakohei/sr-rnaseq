#!/usr/bin/env python3
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Summarize score_normalized per gene x tissue/cancer for GTEx and TCGA")
parser.add_argument("--merged",    required=True, help="merged_results.tsv")
parser.add_argument("--gene_ids",  required=True, help="gene_ids.list (one gene_id per line)")
parser.add_argument("--out_gtex",  required=True, help="Output TSV for GTEx stats")
parser.add_argument("--out_tcga",  required=True, help="Output TSV for TCGA stats")
args = parser.parse_args()

# gene_ids.list の読み込み
with open(args.gene_ids) as f:
    target_genes = set(line.strip() for line in f if line.strip())
print(f"Target genes: {len(target_genes)}")

# merged_results.tsv の読み込み
df = pd.read_csv(args.merged, sep="\t", dtype=str)
print(f"Total rows: {len(df):,}")

# 対象遺伝子に絞る
df = df[df["gene_id"].isin(target_genes)]
print(f"After gene filter: {len(df):,}")

# score_normalized を数値化
df["score_normalized"] = pd.to_numeric(df["score_normalized"], errors="coerce")

# SRP / ERP / DRP で始まらない行（GTEx/TCGAのtissue/cancer projectのみ）
mask_not_sra = ~df["project"].str.match(r"^(SRP|ERP|DRP)", na=False)
df_nonsra = df[mask_not_sra]
print(f"After excluding SRP/ERP/DRP projects: {len(df_nonsra):,}")


def compute_stats(data, source_label, outpath):
    # bw_url が source_label を含む行に絞る
    subset = data[data["bw_url"].str.contains(source_label, case=False, na=False)]
    print(f"  [{source_label}] rows: {len(subset):,} / projects: {subset['project'].nunique()}")

    # 遺伝子 × project ごとに集計
    stats = (
        subset
        .groupby(["gene_id", "project"])["score_normalized"]
        .agg(
            mean   = "mean",
            std    = "std",
            median = "median",
            n      = "count"
        )
        .reset_index()
    )

    stats[["mean", "std", "median"]] = stats[["mean", "std", "median"]].round(2)

    stats.to_csv(outpath, sep="\t", index=False)
    print(f"  Saved to: {outpath}  ({len(stats):,} rows)")


print("\n--- GTEx ---")
compute_stats(df_nonsra, "gtex", args.out_gtex)

print("\n--- TCGA ---")
compute_stats(df_nonsra, "tcga", args.out_tcga)
