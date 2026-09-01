
#!/usr/bin/env python3
"""
Rank genes globally and estimate group counts based on enrichment scores.

This version includes:
- --presence-col : column used for grouping (e.g., Cell type, Cell type group)
- --unique-col   : column used for uniqueness counting
- --unique       : boolean switch to deduplicate based on --unique-col
"""

import argparse
import pandas as pd
import numpy as np


def to_numeric_safe(series: pd.Series, convert_inf_to_nan: bool = True) -> pd.Series:
    out = pd.to_numeric(series, errors="coerce")
    if convert_inf_to_nan:
        out = out.replace([np.inf, -np.inf], np.nan)
    return out


def normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns=lambda x: " ".join(str(x).split()).strip())


def main():
    parser = argparse.ArgumentParser(description="Rank genes globally by enrichment and count groups.")

    # I/O
    parser.add_argument("--input", required=True, help="Input TSV.")
    parser.add_argument("--output", required=True, help="Output TSV.")

    # Top subset selection
    parser.add_argument("--top-percent", type=float, default=25.0)
    parser.add_argument("--min-top-rows", type=int, default=1)
    parser.add_argument("--top-col", default="log2_enrichment_penalized")
    parser.add_argument("--sorting-col", default="log2_enrichment_penalized")

    # Entity columns
    parser.add_argument("--gene-col", default="Gene")
    parser.add_argument("--celltype-col", default="Cell type", help="Deprecated fallback.")
    parser.add_argument("--presence-col", default=None,
                        help="Column defining grouping (e.g., Cell type, Cell type group).")

    # Uniqueness controls
    parser.add_argument("--unique", action="store_true",
                        help="Enable unique counting using --unique-col.")
    parser.add_argument("--unique-col", default=None,
                        help="Column used for uniqueness counting; defaults to --presence-col.")

    # Filters
    parser.add_argument("--drop-na", action="store_true")
    parser.add_argument("--drop-negatives", action="store_true")

    # Output controls
    parser.add_argument("--output-count-col", default="group_count")
    parser.add_argument("--output-list-col", default="group_list")
    parser.add_argument("--output-overall-rank-col", default="overall_rank")
    parser.add_argument("--output-rank-within-col", default="rank_within_group")
    parser.add_argument("--include-cols", nargs="+",
                        default=["Gene", "Gene name", "Cell type", "avg_nCPM",
                                 "specificity_tau", "Enrichment score (tau penalized)",
                                 "log2_enrichment_penalized"])
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    presence_col = args.presence_col if args.presence_col else args.celltype_col
    unique_col = args.unique_col if args.unique_col else presence_col

    # Load file
    df = pd.read_csv(args.input, sep="\t")
    df = normalize_headers(df)

    # Validate
    for col in [args.gene_col, presence_col, unique_col, args.top_col, args.sorting_col]:
        if col not in df.columns:
            raise SystemExit(f"ERROR: Missing column '{col}' in input file.")

    # Numeric processing
    df[args.top_col] = to_numeric_safe(df[args.top_col])
    df[args.sorting_col] = to_numeric_safe(df[args.sorting_col])

    if args.drop_na:
        df = df[df[args.top_col].notna() & df[args.sorting_col].notna()]

    if args.drop_negatives:
        df = df[(df[args.top_col] >= 0) & (df[args.sorting_col] >= 0)]

    if df.empty:
        df.to_csv(args.output, sep="\t", index=False)
        return

    # Global top subset
    df_sorted = df.sort_values(by=args.sorting_col, ascending=False)
    n_top = max(int(len(df_sorted) * (args.top_percent / 100.0)), args.min_top_rows)
    df_top = df_sorted.head(n_top).copy()

    # Counting groups
    grouped = df_top.groupby(args.gene_col)[unique_col]

    if args.unique:
        per_gene_groups = grouped.apply(lambda s: sorted(set(s.dropna())))
    else:
        per_gene_groups = grouped.apply(lambda s: list(s.dropna()))

    per_gene_count = per_gene_groups.apply(len)

    # Attach outputs
    df_top[args.output_count_col] = df_top[args.gene_col].map(per_gene_count).fillna(0).astype("Int64")
    df_top[args.output_list_col] = df_top[args.gene_col].map(
        per_gene_groups.apply(lambda x: ", ".join(x) if isinstance(x, list) else "")
    ).fillna("")

    # Ranking
    df_ranked = df_top.sort_values(
        by=[args.output_count_col, args.sorting_col],
        ascending=[True, False]
    ).reset_index(drop=True)

    df_ranked[args.output_overall_rank_col] = np.arange(1, len(df_ranked) + 1)
    df_ranked[args.output_rank_within_col] = df_ranked.groupby(presence_col).cumcount() + 1

    # Output
    output_cols = args.include_cols + [
        args.output_count_col,
        args.output_list_col,
        args.output_rank_within_col,
        args.output_overall_rank_col,
    ]
    output_cols = [c for c in output_cols if c in df_ranked.columns]

    df_ranked.to_csv(args.output, sep="\t", index=False)

    if args.verbose:
        print(f"[INFO] Wrote {len(df_ranked)} rows to {args.output}")


if __name__ == "__main__":
    main()

