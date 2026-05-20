#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Merge recount3 result files and join with project metadata CSV")
parser.add_argument("--results_dir", required=True, help="Directory containing per-project .txt files")
parser.add_argument("--csv",         required=True, help="recount3 selection CSV (with organism, project, n_samples, study_title, study_abstract columns)")
parser.add_argument("--output",      required=True, help="Output TSV path")
args = parser.parse_args()

# 1. resultsフォルダのtxtを全結合（project列を追加）
print("Loading result files...")
dfs = []
for path in glob.glob(os.path.join(args.results_dir, "*.txt")):
    project = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "project", project)
    dfs.append(df)

combined = pd.concat(dfs, ignore_index=True)
print(f"  Combined: {len(combined):,} rows from {len(dfs)} files")

# 2. CSVからhuman行だけ抽出し、必要列を選択
meta = pd.read_csv(args.csv)
human_meta = meta[meta["organism"] == "human"][["project", "n_samples", "study_title", "study_abstract"]]
print(f"  Human projects in CSV: {len(human_meta)}")

# 3. left join（combinedの行は削除しない）
result = combined.merge(human_meta, on="project", how="left")
print(f"  After merge: {len(result):,} rows")

out_dir = os.path.dirname(args.output)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

result.to_csv(args.output, sep="\t", index=False)
print(f"Saved to: {args.output}")
