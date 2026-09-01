
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract tokens from a delimited column that match a regex pattern.

Features:
- Regex-based matching (prefix by default via re.match, or anywhere via re.search).
- Extract the nth match (supports multiple nth values, e.g., --nth 1 2 3).
- Or extract ALL matches in a single column via --all-matches.
- Handles large files with --chunksize streaming.
- Fully configurable delimiters and encoding.

Examples:
    # First value that starts with Predicted or Secreted
    python extract_regex.py \
        --input genes.tsv \
        --output out.tsv \
        --file-delim $'\t' \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern 'Predicted|Secreted' \
        --ignore-case \
        --nth 1

    # First and second matches + all matches
    python extract_regex.py \
        --input genes.tsv \
        --output out_all.tsv \
        --file-delim $'\t' \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern 'Predicted|Secreted' \
        --ignore-case \
        --nth 1 2 \
        --all-matches \
        --join-matches-sep "; "
"""

import argparse
import sys
import re
from typing import List, Optional
import pandas as pd


# ---------------------------
# Core helpers
# ---------------------------

def split_tokens(cell: Optional[str], sep: str) -> List[str]:
    """Split a cell string by sep and trim whitespace. Return [] for null/empty."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    text = str(cell).strip()
    if not text:
        return []
    return [tok.strip() for tok in text.split(sep)]


def compile_regex(pattern: str, ignore_case: bool) -> re.Pattern:
    """Compile the regex with optional IGNORECASE; raise ValueError on invalid regex."""
    flags = re.IGNORECASE if ignore_case else 0
    try:
        return re.compile(pattern, flags)
    except re.error as rex:
        raise ValueError(f"Invalid regex pattern '{pattern}': {rex}")


def find_matches_in_tokens(tokens: List[str], rx: re.Pattern, match_anywhere: bool) -> List[str]:
    """Return tokens that match the regex (match at start or anywhere depending on flag)."""
    matcher = rx.search if match_anywhere else rx.match
    return [t for t in tokens if matcher(t)]


# ---------------------------
# Processing
# ---------------------------

def process_chunk(df: pd.DataFrame,
                  id_cols: List[str],
                  target_col: str,
                  pattern: str,
                  value_sep: str,
                  nth_list: Optional[List[int]],
                  all_matches: bool,
                  join_matches_sep: str,
                  ignore_case: bool,
                  match_anywhere: bool) -> pd.DataFrame:
    """Process a DataFrame chunk and return the output DataFrame with requested columns."""
    # Validate columns
    if target_col not in df.columns:
        raise KeyError(f"Target column '{target_col}' not found. Available: {list(df.columns)}")

    missing_id = [c for c in id_cols if c not in df.columns]
    if missing_id:
        raise KeyError(f"ID columns not found: {missing_id}. Available: {list(df.columns)}")

    # Determine behavior: default to nth=[1] if neither nth nor all-matches provided
    if (not all_matches) and (not nth_list):
        nth_list = [1]

    # Validate nth_list (positive integers) and de-duplicate while preserving order
    if nth_list:
        seen = set()
        clean_nth = []
        for n in nth_list:
            if not isinstance(n, int) or n <= 0:
                raise ValueError(f"--nth values must be positive integers. Got: {n}")
            if n not in seen:
                clean_nth.append(n)
                seen.add(n)
        nth_list = clean_nth

    # Compile regex once per chunk
    rx = compile_regex(pattern, ignore_case)

    # Compute matches once per row
    matches_series = df[target_col].apply(
        lambda cell: find_matches_in_tokens(split_tokens(cell, value_sep), rx, match_anywhere)
    )

    # Start output with ID columns
    out_df = df[id_cols].copy()

    # Add nth columns
    if nth_list:
        for n in nth_list:
            colname = f"{target_col}__regex_nth{n}"
            out_df[colname] = matches_series.apply(lambda m: m[n - 1] if len(m) >= n else "")

    # Add all-matches column
    if all_matches:
        colname_all = f"{target_col}__regex_all"
        out_df[colname_all] = matches_series.apply(lambda m: join_matches_sep.join(m) if m else "")

    return out_df


# ---------------------------
# CLI
# ---------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract nth or all regex-matching tokens from a delimited column."
    )
    p.add_argument("--input", required=True, help="Path to input file (TSV/CSV/other delimited).")
    p.add_argument("--output", required=True, help="Path to output file.")
    p.add_argument("--file-delim", default="\t", help=r"Input file delimiter. Default: '\t'.")
    p.add_argument("--output-delim", default="\t", help=r"Output file delimiter. Default: '\t'.")
    p.add_argument("--encoding", default="utf-8", help="Input file encoding. Default: utf-8.")
    p.add_argument("--id-cols", nargs="+", required=True,
                   help="Columns to include in output (e.g., Gene Ensembl).")
    p.add_argument("--target-col", required=True,
                   help="Column containing delimited values (e.g., 'Protein class').")
    p.add_argument("--value-sep", default=",",
                   help="Separator used inside the target column values. Default: ','.")
    p.add_argument("--pattern", required=True,
                   help=("Regex pattern to match tokens. Examples: 'Predicted|Secreted' (prefix by default), "
                         "'^(Predicted|Secreted)', or '.*(Predicted|Secreted).*' with --match-anywhere."))
    p.add_argument("--nth", nargs="+", type=int, default=None,
                   help="One or more 1-based match indexes to extract (e.g., --nth 1 2 3).")
    p.add_argument("--all-matches", action="store_true",
                   help="Also output ALL matching tokens joined into one column.")
    p.add_argument("--join-matches-sep", default="; ",
                   help="Separator used when joining all matches. Default: '; '.")
    p.add_argument("--ignore-case", action="store_true",
                   help="Case-insensitive regex matching.")
    p.add_argument("--match-anywhere", action="store_true",
                   help="Use re.search() (match anywhere). Default uses re.match() (prefix only).")
    p.add_argument("--chunksize", type=int, default=None,
                   help="Read/process file in chunks (rows per chunk) for large datasets.")
    return p


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    try:
        if args.chunksize:
            wrote_header = False
            for chunk in pd.read_csv(
                args.input,
                sep=args.file_delim,
                dtype=str,
                chunksize=args.chunksize,
                encoding=args.encoding,
                keep_default_na=False
            ):
                out_chunk = process_chunk(
                    chunk,
                    id_cols=args.id_cols,
                    target_col=args.target_col,
                    pattern=args.pattern,
                    value_sep=args.value_sep,
                    nth_list=args.nth,
                    all_matches=args.all_matches,
                    join_matches_sep=args.join_matches_sep,
                    ignore_case=args.ignore_case,
                    match_anywhere=args.match_anywhere
                )
                out_chunk.to_csv(
                    args.output,
                    sep=args.output_delim,
                    index=False,
                    mode="w" if not wrote_header else "a",
                    header=not wrote_header,
                    encoding="utf-8"
                )
                wrote_header = True
        else:
            df = pd.read_csv(
                args.input,
                sep=args.file_delim,
                dtype=str,
                encoding=args.encoding,
                keep_default_na=False
            )
            out_df = process_chunk(
                df,
                id_cols=args.id_cols,
                target_col=args.target_col,
                pattern=args.pattern,
                value_sep=args.value_sep,
                nth_list=args.nth,
                all_matches=args.all_matches,
                join_matches_sep=args.join_matches_sep,
                ignore_case=args.ignore_case,
                match_anywhere=args.match_anywhere
            )
            out_df.to_csv(args.output, sep=args.output_delim, index=False, encoding="utf-8")

    except KeyError as ke:
        print(f"[Error] {ke}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as fe:
        print(f"[Error] {fe}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.ParserError as pe:
        print(f"[Parse Error] {pe}", file=sys.stderr)
        sys.exit(1)
    except ValueError as ve:
        print(f"[Regex/Argument Error] {ve}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[Unexpected Error] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()