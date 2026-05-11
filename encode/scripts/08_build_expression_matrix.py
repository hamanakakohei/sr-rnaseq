"""
Build expression matrix:
- Rows: RNA-seq experiments (ENCSR) with metadata
- Columns: gene_id x strand, ordered as gene1_plus, gene1_minus, gene2_plus, gene2_minus, ...
- Values: total_signal from bigWig expression files
"""
import argparse
import pandas as pd
from pathlib import Path

META_COLS = [
    "Biosample classification",
    "Biosample summary",
    "Biosample term name",
    "Dbxrefs",
    "Description",
    "Organism",
    "Life stage",
    "Biosample age",
    "Cellular component",
    "Biosample ontology",
    "Simple Biosample summary",
    "Life stage and age summary",
]

def load_expression(expr_dir, encff_id):
    path = expr_dir / f"{encff_id}.TotalSignal.txt"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return df.set_index("gene_id")["total_signal"]

def main():
    parser = argparse.ArgumentParser(description="Build a gene_id x ENCSR expression matrix with metadata.")
    parser.add_argument("--expr-dir", type=Path, help="Directory containing ENCFF*.TotalSignal.txt files (default: %(default)s)")
    parser.add_argument("--mapping", type=Path, help="TSV mapping ENCSR to plus/minus ENCFF IDs (default: %(default)s)")
    parser.add_argument("--summary", type=Path, help="TSV with ENCSR metadata (default: %(default)s)")
    parser.add_argument("--output", type=Path, help="Output TSV path (default: %(default)s)")
    args = parser.parse_args()

    # Load ENCSR -> (minus ENCFF, plus ENCFF) mapping
    mapping = pd.read_csv(args.mapping, sep="\t")
    print(f"ENCsrs: {len(mapping)}")

    # Load summary metadata keyed by encsr (deduplicate)
    summary = pd.read_csv(args.summary, sep="\t", low_memory=False)
    meta_cols_present = [c for c in META_COLS if c in summary.columns]
    summary_meta = (
        summary[["encsr"] + meta_cols_present]
        .drop_duplicates(subset="encsr")
        .set_index("encsr")
    )
    print(f"ENCsrs with metadata: {len(summary_meta)}")

    # Determine gene_id order from first available file
    first_encff = mapping["plus"].iloc[0]
    first_series = load_expression(args.expr_dir, first_encff)
    if first_series is None:
        raise FileNotFoundError(f"Could not load {first_encff}")
    gene_ids = list(first_series.index)
    print(f"Gene IDs: {len(gene_ids)}")

    # Build interleaved column order: gene1_plus, gene1_minus, gene2_plus, gene2_minus, ...
    col_order = []
    for g in gene_ids:
        col_order.append(f"{g}_plus")
        col_order.append(f"{g}_minus")

    rows = []
    for _, row in mapping.iterrows():
        encsr = row["encsr"]
        plus_encff = row["plus"]
        minus_encff = row["minus"]

        plus_expr = load_expression(args.expr_dir, plus_encff)
        minus_expr = load_expression(args.expr_dir, minus_encff)

        record = {"encsr": encsr}

        # Add metadata
        if encsr in summary_meta.index:
            for col in meta_cols_present:
                record[col] = summary_meta.at[encsr, col]
        else:
            for col in meta_cols_present:
                record[col] = pd.NA

        # Add expression values
        for g in gene_ids:
            record[f"{g}_plus"] = plus_expr[g] if (plus_expr is not None and g in plus_expr.index) else pd.NA
            record[f"{g}_minus"] = minus_expr[g] if (minus_expr is not None and g in minus_expr.index) else pd.NA

        rows.append(record)
        if len(rows) % 20 == 0:
            print(f"  Processed {len(rows)}/{len(mapping)} experiments...")

    df = pd.DataFrame(rows).set_index("encsr")

    # Reorder: metadata first, then interleaved expression columns
    meta_cols_in_df = [c for c in meta_cols_present if c in df.columns]
    df = df[meta_cols_in_df + col_order]

    df.to_csv(args.output, sep="\t")
    print(f"\nDone. Output: {args.output}")
    print(f"Shape: {df.shape} (rows x cols)")

if __name__ == "__main__":
    main()
