
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Match combinations of columns between two tables, collect values from the second table,
and attach them to the first table as a new column.

Example:
    python match_and_collect.py \
      file1.xlsx file2.tsv \
      --file1-keys "Gene name" "Cell type" \
      --file2-keys "Gene name" "Cell type" \
      --collect-column "Tissue" \
      --filter-col nCPM --filter-op ">" --filter-value 0 \
      --output-column "Matched tissues" \
      --sep-out "; " \
      --case-insensitive \
      --drop-dupes \
      --output merged.xlsx
"""

import argparse
import os
import sys
from typing import List, Optional
import pandas as pd


def read_table(path: str, sep: Optional[str] = None, sheet: Optional[str] = None) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in [".csv"]:
        return pd.read_csv(path, sep=sep or ",")
    elif ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep=sep or "\t")
    elif ext in [".xlsx", ".xls"]:
        # Prefer explicit engines for compatibility
        if ext == ".xlsx":
            try:
                return pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
            except Exception:
                # Fallback if engine is unavailable in the environment
                return pd.read_excel(path, sheet_name=sheet)
        else:
            try:
                return pd.read_excel(path, sheet_name=sheet, engine="xlrd")
            except Exception:
                return pd.read_excel(path, sheet_name=sheet)
    else:
        raise ValueError(f"Unsupported file extension for {path!r}. "
                         f"Use .csv, .tsv, .txt, .xlsx, or .xls")


def coerce_numeric_if_possible(series: pd.Series):
    """Try convert to numeric while leaving non-convertible as NaN."""
    return pd.to_numeric(series, errors="coerce")


def apply_simple_filter(df: pd.DataFrame, col: str, op: str, val: str) -> pd.DataFrame:
    if col not in df.columns:
        raise KeyError(f"Filter column {col!r} not found in the second file.")

    s = df[col]
    # Try to interpret the value as numeric; if conversion fails, treat as string
    try:
        v_num = float(val)
        s_num = coerce_numeric_if_possible(s)
        if s_num.isna().all():
            # Column couldn't be numeric; fall back to string comparison
            raise ValueError
        if op == ">":
            mask = s_num > v_num
        elif op == ">=":
            mask = s_num >= v_num
        elif op == "<":
            mask = s_num < v_num
        elif op == "<=":
            mask = s_num <= v_num
        elif op == "==":
            mask = s_num == v_num
        elif op == "!=":
            mask = s_num != v_num
        else:
            raise ValueError(f"Unsupported operator {op!r}. Use one of: > >= < <= == !=")
        return df[mask]
    except Exception:
        # String comparison path
        if op not in ["==", "!="]:
            raise ValueError(f"Operator {op!r} is only valid for numeric comparisons unless you use == or !=.")
        if op == "==":
            mask = s.astype(str) == str(val)
        else:
            mask = s.astype(str) != str(val)
        return df[mask]


def normalize_keys(df: pd.DataFrame, key_cols: List[str], case_insensitive: bool, strip_ws: bool) -> pd.DataFrame:
    df = df.copy()
    for c in key_cols:
        if c not in df.columns:
            raise KeyError(f"Key column {c!r} not found.")
        col = df[c].astype(str)
        if strip_ws:
            col = col.str.strip()
        if case_insensitive:
            col = col.str.lower()
        df[c] = col
    return df


def aggregate_values(
    df2: pd.DataFrame,
    file2_keys: List[str],
    collect_column: str,
    dedupe: bool,
    sort_values: bool,
    sep_out: str,
    output_colname: str
) -> pd.DataFrame:
    if collect_column not in df2.columns:
        raise KeyError(f"Collect column {collect_column!r} not found in the second file.")
    # Keep only the keys + collect column to avoid suffix collisions later
    sub = df2[file2_keys + [collect_column]].copy()
    # Group and aggregate
    def agg_func(s: pd.Series) -> str:
        vals = s.dropna().astype(str).tolist()
        if dedupe:
            # Preserve first-seen order when deduping unless sorting requested
            seen = set()
            vals = [x for x in vals if (x not in seen and not seen.add(x))]
        if sort_values:
            vals = sorted(vals)
        return sep_out.join(vals)

    grouped = (
        sub.groupby(file2_keys, dropna=False)[collect_column]
           .apply(agg_func)
           .reset_index()
           .rename(columns={collect_column: output_colname})
    )
    return grouped


def main():
    parser = argparse.ArgumentParser(
        description="Match column combinations between two files, collect values from the second, and attach to the first."
    )
    parser.add_argument("file1", help="Path to the first file (CSV/TSV/TXT/XLSX/XLS).")
    parser.add_argument("file2", help="Path to the second file (CSV/TSV/TXT/XLSX/XLS).")
    parser.add_argument("--file1-sep", help="Delimiter for file1 if CSV/TXT/TSV (default: auto by extension).")
    parser.add_argument("--file2-sep", help="Delimiter for file2 if CSV/TXT/TSV (default: auto by extension).")
    parser.add_argument("--sheet1", help="Excel sheet name for file1 (if .xlsx/.xls).")
    parser.add_argument("--sheet2", help="Excel sheet name for file2 (if .xlsx/.xls).")

    parser.add_argument("--file1-keys", nargs="+", required=True, help="Key columns in file1 (e.g., 'Gene name' 'Cell type').")
    parser.add_argument("--file2-keys", nargs="+", required=True, help="Corresponding key columns in file2 (same count and order).")
    parser.add_argument("--collect-column", required=True, help="Column from file2 to collect (e.g., 'Tissue').")

    # Filtering options (either a simple triplet OR a pandas.query expression)
    parser.add_argument("--filter-col", help="Filter column in file2 (e.g., nCPM).")
    parser.add_argument("--filter-op", help="Filter operator: >, >=, <, <=, ==, !=")
    parser.add_argument("--filter-value", help="Filter value (e.g., 0).")
    parser.add_argument("--filter-expr", help="Advanced pandas.query expression on file2 (e.g., \"`nCPM` > 0 and `Cell count` >= 10\").")

    parser.add_argument("--output-column", required=True, help="Name of the new column to add to file1 (e.g., 'Matched tissues').")
    parser.add_argument("--sep-out", default="; ", help="Separator for collected values (default: '; ').")
    parser.add_argument("--drop-dupes", action="store_true", help="De-duplicate collected values.")
    parser.add_argument("--sort-values", action="store_true", help="Sort collected values alphabetically.")
    parser.add_argument("--case-insensitive", action="store_true", help="Lowercase key values before matching.")
    parser.add_argument("--no-strip", action="store_true", help="Disable trimming whitespace in key columns.")
    parser.add_argument("--output", help="Where to save the merged result (CSV/XLSX). If omitted, prints a preview.")

    args = parser.parse_args()

    if len(args.file1_keys) != len(args.file2_keys):
        print("Error: --file1-keys and --file2-keys must have the same number of columns (same order).", file=sys.stderr)
        sys.exit(1)

    # Read files
    try:
        df1 = read_table(args.file1, sep=args.file1_sep, sheet=args.sheet1)
        df2 = read_table(args.file2, sep=args.file2_sep, sheet=args.sheet2)
    except Exception as e:
        print(f"Failed to read input files: {e}", file=sys.stderr)
        sys.exit(1)

    strip_ws = not args.no_strip

    # Normalize keys
    try:
        df1 = normalize_keys(df1, args.file1_keys, case_insensitive=args.case_insensitive, strip_ws=strip_ws)
        df2 = normalize_keys(df2, args.file2_keys, case_insensitive=args.case_insensitive, strip_ws=strip_ws)
    except KeyError as e:
        print(f"Key error: {e}", file=sys.stderr)
        sys.exit(1)

    # Optional filtering on df2
    if args.filter_expr and any([args.filter_col, args.filter_op, args.filter_value]):
        print("Error: Use either --filter-expr OR the simple --filter-col/--filter-op/--filter-value triplet, not both.", file=sys.stderr)
        sys.exit(1)

    if args.filter_expr:
        try:
            # Allow backticks around columns with spaces
            df2 = df2.query(args.filter_expr, engine="python")
        except Exception as e:
            print(f"Failed to apply filter expression: {e}", file=sys.stderr)
            sys.exit(1)
    elif all([args.filter_col, args.filter_op, args.filter_value is not None]):
        try:
            df2 = apply_simple_filter(df2, args.filter_col, args.filter_op, args.filter_value)
        except Exception as e:
            print(f"Failed to apply simple filter: {e}", file=sys.stderr)
            sys.exit(1)
    elif any([args.filter_col, args.filter_op, args.filter_value is not None]):
        print("Error: Incomplete filter specification. Provide all of --filter-col, --filter-op, and --filter-value; or use --filter-expr.", file=sys.stderr)
        sys.exit(1)

    # Aggregate values in df2 by keys
    try:
        agg = aggregate_values(
            df2=df2,
            file2_keys=args.file2_keys,
            collect_column=args.collect_column,
            dedupe=args.drop_dupes,
            sort_values=args.sort_values,
            sep_out=args.sep_out,
            output_colname=args.output_column
        )
    except Exception as e:
        print(f"Aggregation failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Merge onto df1
    try:
        merged = df1.merge(
            agg,
            how="left",
            left_on=args.file1_keys,
            right_on=args.file2_keys
        )
    except Exception as e:
        print(f"Merge failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Remove duplicate key columns from right side if they differ in name
    # (not necessary here since 'agg' only has keys + new column)
    # Ensure the result keeps original file1 columns + new output column at the end
    # Already the case with the merge above.

    if args.output:
        out_ext = os.path.splitext(args.output)[1].lower()
        try:
            if out_ext in [".xlsx", ".xls"]:
                merged.to_excel(args.output, index=False, engine=("openpyxl" if out_ext == ".xlsx" else None))
            else:
                # Default to CSV; allow user to specify a custom delimiter by file name if they want .tsv
                sep = "\t" if out_ext in [".tsv", ".txt"] else ","
                merged.to_csv(args.output, index=False, sep=sep)
            print(f"Saved merged table to: {args.output}")
        except Exception as e:
            print(f"Failed to save output: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Print a concise preview
        print("\nMerged preview (top 10 rows):")
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(merged.head(10))
        print("\n(No --output provided; nothing was written to disk.)")


if __name__ == "__main__":
    main()