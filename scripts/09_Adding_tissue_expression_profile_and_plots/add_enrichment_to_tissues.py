#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Propagate enrichment values from table1 into tissue strings in table2."
    )

    # I/O
    parser.add_argument("--table1", required=True, help="TSV with Tissue and enrichment values")
    parser.add_argument("--table2", required=True, help="TSV with Present tissues column")
    parser.add_argument("--output", required=True, help="Output TSV")

    # Join keys
    parser.add_argument("--gene-col", default="Gene")
    parser.add_argument("--gene-name-col", default="Gene name")
    parser.add_argument("--cell-type-col", default="Cell type")

    # Tissue + enrichment columns
    parser.add_argument("--tissue-col-table1", default="Tissue",
                        help="Column in table1 containing single tissue name")
    parser.add_argument("--tissue-col-table2", default="Present tissues",
                        help="Column in table2 containing delimited tissue names")
    parser.add_argument("--enrichment-col", default="log2_enrichment_penalized",
                        help="Column in table1 with the value to propagate")

    # Formatting / behavior
    parser.add_argument("--tissue-sep", default="&",
                        help="Delimiter used between tissues in table2 (default: '&')")
    parser.add_argument("--value-sep", default=":",
                        help="Separator between value and tissue name (default: ':')")
    parser.add_argument("--missing-value", default="NA",
                        help="Value/string to use when a tissue has no match in table1")
    parser.add_argument("--case-insensitive", action="store_true",
                        help="Case-insensitive matching for keys (gene, gene name, cell type, tissue)")
    parser.add_argument("--output-tissues-col", default=None,
                        help="If provided, write the enriched tissues to this NEW column instead of "
                             "overwriting the original tissue column")

    # NEW: preprocess tissue from table1 for matching (e.g., 'breast-1' -> 'breast')
    parser.add_argument("--tissue1-trim-after", default=None,
                        help="Split tissue at the FIRST occurrence of this string and keep the part BEFORE it. "
                             "Example: --tissue1-trim-after '-' turns 'breast-1' into 'breast'.")
    parser.add_argument("--tissue1-regex", default=None,
                        help="Regex to EXTRACT the portion of the table1 tissue to match. "
                             "Example: --tissue1-regex '^[^-]+' turns 'breast-1' into 'breast'. "
                             "Takes precedence over --tissue1-trim-after.")

    # NEW: choose which tissue label to render in the output
    parser.add_argument("--label-source", choices=["table2", "table1"], default="table2",
                        help="Which tissue label to print in the enriched output: "
                             "'table2' (default) prints the token from table2; "
                             "'table1' prints the ORIGINAL table1 tissue when matched, "
                             "falling back to table2 token if unmatched.")

    return parser.parse_args()


def norm(x, lower=False):
    """Normalize scalars for matching: ensure string, strip whitespace, optionally lowercase."""
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return s.lower() if lower else s


def preprocess_tissue1(raw_tissue: str, lower: bool, trim_after: str | None, regex: str | None) -> str:
    """Apply optional preprocessing to Table-1 tissue value before using as a key."""
    s = norm(raw_tissue, lower=False)  # normalize whitespace but do NOT lowercase yet
    if regex:
        m = re.search(regex, s)
        s = m.group(0) if m else s
    elif trim_after:
        idx = s.find(trim_after)
        if idx != -1:
            s = s[:idx]
    # final normalization and optional case-fold
    return norm(s, lower)


def main():
    args = parse_args()

    # Load as strings; do not let pandas guess types
    try:
        t1 = pd.read_csv(args.table1, sep="\t", dtype=str, keep_default_na=False)
        t2 = pd.read_csv(args.table2, sep="\t", dtype=str, keep_default_na=False)
    except Exception as e:
        print(f"[ERROR] Failed to read input TSV(s): {e}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns exist
    needed_t1 = {args.gene_col, args.gene_name_col, args.cell_type_col,
                 args.tissue_col_table1, args.enrichment_col}
    missing_t1 = [c for c in needed_t1 if c not in t1.columns]
    if missing_t1:
        print(f"[ERROR] Missing columns in table1: {missing_t1}", file=sys.stderr)
        sys.exit(1)

    needed_t2 = {args.gene_col, args.gene_name_col, args.cell_type_col, args.tissue_col_table2}
    missing_t2 = [c for c in needed_t2 if c not in t2.columns]
    if missing_t2:
        print(f"[ERROR] Missing columns in table2: {missing_t2}", file=sys.stderr)
        sys.exit(1)

    # Whether to match case-insensitive
    lower = args.case_insensitive

    # Build lookup dict from table1
    # Key: (gene, gene_name, cell_type, tissue_processed) -> dict(value, t1_label)
    lookup = {}
    for _, row in t1.iterrows():
        t1_raw_label = norm(row[args.tissue_col_table1], lower=False)  # keep original label
        tissue_processed = preprocess_tissue1(
            t1_raw_label,
            lower=lower,
            trim_after=args.tissue1_trim_after,
            regex=args.tissue1_regex,
        )
        key = (
            norm(row[args.gene_col], lower),
            norm(row[args.gene_name_col], lower),
            norm(row[args.cell_type_col], lower),
            tissue_processed,
        )
        lookup[key] = {
            "value": norm(row[args.enrichment_col], lower=False),
            "t1_label": t1_raw_label,
        }
        # Note: if duplicates map to the same key, the last one wins.
        # Sort table1 beforehand if you need deterministic precedence.

    # Prepare output column name
    out_col = args.tissue_col_table2 if args.output_tissues_col is None else args.output_tissues_col

    # Enrichment function for each row in table2
    def enrich_row(row):
        gene = norm(row[args.gene_col], lower)
        gname = norm(row[args.gene_name_col], lower)
        ctype = norm(row[args.cell_type_col], lower)

        raw_tissues = row[args.tissue_col_table2]
        if pd.isna(raw_tissues) or raw_tissues == "":
            return raw_tissues  # nothing to do

        # Split by the provided separator and trim spaces around each token
        tissues = [t.strip() for t in str(raw_tissues).split(args.tissue_sep)]
        enriched_parts = []

        for t2_token in tissues:
            t2_key = norm(t2_token, lower)  # normalized table2 token for matching
            key = (gene, gname, ctype, t2_key)

            if key in lookup:
                value = lookup[key]["value"]
                label = lookup[key]["t1_label"] if args.label_source == "table1" else t2_token
            else:
                value = args.missing_value
                label = t2_token  # fall back to table2 token when no match

            enriched_parts.append(f"{value}{args.value_sep}{label}")

        # Re-join with the same spacing convention: " <sep> "
        return f" {args.tissue_sep} ".join(enriched_parts)

    # Compute enriched tissues
    t2[out_col] = t2.apply(enrich_row, axis=1)

    # Write output
    try:
        t2.to_csv(args.output, sep="\t", index=False)
    except Exception as e:
        print(f"[ERROR] Failed to write output TSV: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
