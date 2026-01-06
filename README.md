# cell_type_enrichment_v2
This is an enhanced version of cell type enrichment script to identify cell type specific gene expression and cell type enriched gene expression

This pipeline uses publically available data.

I use the datasets from Human protein Atlas as they have a extensive dataset they create by combinitng multiple studies

You can find them in
```url
https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data
```

#*First of all, you should install needed python packages and load them*

I have done this and I load needed packages in alliancecan like following
```bash
module load gcc arrow
module load python
python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate
```

# 1) Combining datasets

I am using two main datasets available here

HPA has,

•  Raw counts → pseudobulk within each cluster
	They average raw counts across cells inside each cluster, producing one pseudobulk count vector per cluster.
	
•  pCPM: counts per million normalization
	CPM is computed from pseudobulk counts (not per cell counts).
	CPM does not include a weighting by number of cells — it is based only on total pseudobulk counts.
	
•  TMM normalization of pCPM → nCPM
	TMM normalization scales CPM profiles to remove compositional bias.
	TMM does not use cell counts; it uses global expression distributions.

These nCPM values are in “rna_single_cell_cluster.tsv”

```tsv
Gene    Gene name       Tissue  Cluster Cell type       Read count      nCPM
ENSG00000000003 TSPAN6  ovary   c-0     ovarian stromal cells   493     92.5
ENSG00000000003 TSPAN6  ovary   c-1     ovarian stromal cells   529     80.5
ENSG00000000003 TSPAN6  ovary   c-2     ovarian stromal cells   143     52.3
ENSG00000000003 TSPAN6  ovary   c-3     ovarian stromal cells   456     91.4
ENSG00000000003 TSPAN6  ovary   c-4     vascular endothelial cells      164     28.6
```

The file “rna_single_cell_clusters.tsv” (which is different from the file above) contains information on cell counts and reliability
```tsv
Tissue  Cluster Cell type       Cell type detail        Cell type class Cell count      Included in aggregation Annotation reliability
adipose tissue  c-0     mesothelial cells       mesothelial cells       specialized epithelial cells    8942    yes     high
adipose tissue  c-1     adipocytes      mature adipocytes       mesenchymal cells       6996    yes     high
adipose tissue  c-2     adipocytes      mature adipocytes       mesenchymal cells       6993    yes     high
```
Because I need cell count data for my analysis, I start by combining these dataframes to add columns "Cell count" , "Included in aggregation", "Annotation reliability" from rna_single_cell_clusters.tsv to rna_single_cell_cluster.tsv matching by Cluster.
I did this with merge_tsv_by_keys.py

merge_tsv_by_keys.py
```py

#!/usr/bin/env python3
"""
Merge any two TSV files by user-specified key columns and append selected columns
from the right file into the left file.

Features:
- Arbitrary merge keys: same names or mapped across files.
- Default separator is tab; configurable via CLI.
- Select which columns to append from the right file (default: all non-key columns).
- Left join by default; supports inner/right/outer.
- Optional merge validation (m:1, 1:1, m:m).

Usage examples:
    python merge_tsv_by_keys.py \
        --left rna_single_cell_cluster.tsv \
        --right rna_single_cell_clusters.tsv \
        --left-keys Cluster \
        --right-keys Cluster \
        --right-cols "Cell count,Included in aggregation,Annotation reliability" \
        --out enriched.tsv

    # When key names differ between files:
    python merge_tsv_by_keys.py \
        --left left.tsv --right right.tsv \
        --map "cluster_id:Cluster,tissue_name:Tissue" \
        --out merged.tsv
"""

import argparse
import sys
import os
import pandas as pd


def parse_comma_list(s: str):
    return [x.strip() for x in s.split(",")] if s else []


def parse_key_mapping(s: str):
    """
    Parse a mapping string like: "leftKey1:rightKeyA,leftKey2:rightKeyB"
    Returns (left_keys, right_keys) lists in aligned order.
    """
    if not s:
        return [], []
    left_keys, right_keys = [], []
    for pair in s.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if ":" not in pair:
            raise ValueError(f"Invalid mapping entry '{pair}'. Expected 'leftKey:rightKey'.")
        l, r = [p.strip() for p in pair.split(":", 1)]
        if not l or not r:
            raise ValueError(f"Invalid mapping entry '{pair}'. Empty side.")
        left_keys.append(l)
        right_keys.append(r)
    return left_keys, right_keys


def main():
    parser = argparse.ArgumentParser(
        description="Merge two TSVs on arbitrary key columns and append selected columns from the right file."
    )
    parser.add_argument("--left", required=True, help="Filename of the left/base TSV (will be enriched).")
    parser.add_argument("--right", required=True, help="Filename of the right/metadata TSV (columns to append).")
    parser.add_argument("--out", default=None, help="Output filename (TSV). Default: <left_basename>_merged.tsv")

    # Key specification (two ways):
    parser.add_argument("--left-keys", default=None,
                        help="Comma-separated key columns in LEFT file (e.g., 'Cluster' or 'Tissue,Cluster').")
    parser.add_argument("--right-keys", default=None,
                        help="Comma-separated key columns in RIGHT file (must align with --left-keys).")
    parser.add_argument("--map", default=None,
                        help="Alternative to --left-keys/--right-keys. "
                             "Map key names across files: 'leftKey:rightKey,leftKey2:rightKey2'.")

    # Columns to append from RIGHT (default: all non-key columns)
    parser.add_argument("--right-cols", default=None,
                        help="Comma-separated columns from RIGHT to append. "
                             "Default: all non-key columns in RIGHT.")

    # Merge behavior
    parser.add_argument("--how", choices=["left", "inner", "right", "outer"], default="left",
                        help="Join type (default: left).")
    parser.add_argument("--validate", choices=["m:1", "1:1", "m:m"], default="m:1",
                        help="Merge validation (default: m:1).")
    parser.add_argument("--sep", default="\t",
                        help="Field separator for input/output. Default: tab ('\\t').")

    args = parser.parse_args()

    left_file = os.path.basename(args.left)
    right_file = os.path.basename(args.right)
    out_file = args.out or f"{os.path.splitext(left_file)[0]}_merged.tsv"

    # Determine keys
    if args.map:
        left_keys, right_keys = parse_key_mapping(args.map)
        if not left_keys or not right_keys:
            sys.exit("ERROR: --map provided but no valid key pairs found.")
    else:
        if not args.left_keys or not args.right_keys:
            sys.exit("ERROR: Specify keys via --left-keys and --right-keys, or use --map.")
        left_keys = parse_comma_list(args.left_keys)
        right_keys = parse_comma_list(args.right_keys)
        if len(left_keys) != len(right_keys):
            sys.exit("ERROR: --left-keys and --right-keys must have the same number of columns.")

    # Read files
    try:
        df_left = pd.read_csv(left_file, sep=args.sep, dtype=str)
    except Exception as e:
        sys.exit(f"Error reading LEFT file '{left_file}': {e}")
    try:
        df_right = pd.read_csv(right_file, sep=args.sep, dtype=str)
    except Exception as e:
        sys.exit(f"Error reading RIGHT file '{right_file}': {e}")

    # Validate key columns exist
    missing_left = [c for c in left_keys if c not in df_left.columns]
    missing_right = [c for c in right_keys if c not in df_right.columns]
    if missing_left:
        sys.exit(f"ERROR: LEFT file '{left_file}' missing key columns: {missing_left}")
    if missing_right:
        sys.exit(f"ERROR: RIGHT file '{right_file}' missing key columns: {missing_right}")

    # Decide which right columns to append
    if args.right_cols:
        right_cols_to_add = parse_comma_list(args.right_cols)
        missing_rc = [c for c in right_cols_to_add if c not in df_right.columns]
        if missing_rc:
            sys.exit(f"ERROR: RIGHT file '{right_file}' missing requested columns: {missing_rc}")
    else:
        # Default: all non-key columns from RIGHT
        right_cols_to_add = [c for c in df_right.columns if c not in right_keys]
        if not right_cols_to_add:
            sys.exit("ERROR: No non-key columns found in RIGHT to append. Use --right-cols to specify.")

    # Trim RIGHT to keys + selected columns and drop duplicate key rows
    df_right_trim = df_right[right_keys + right_cols_to_add].drop_duplicates(subset=right_keys)

    # Perform merge
    try:
        df_out = df_left.merge(
            df_right_trim,
            left_on=left_keys,
            right_on=right_keys,
            how=args.how,
            validate=args.validate
        )
    except Exception as e:
        sys.exit(f"Merge error: {e}")

    # If key names differ, pandas will keep both sets; we can optionally drop the right key duplicates.
    # Here we drop the right-side key columns if they duplicate left-side names.
    # If the names differ, both will be kept (useful for auditing).
    for lk, rk in zip(left_keys, right_keys):
        if lk == rk:
            # Same name — pandas keeps one, so nothing to drop
            continue
        # Different names — keep both for transparency (comment below to drop)
        # If you prefer to drop RIGHT keys after merge, uncomment:
        # df_out.drop(columns=[rk], inplace=True, errors="ignore")

    # Write output
    try:
        df_out.to_csv(out_file, sep=args.sep, index=False)
    except Exception as e:
        sys.exit(f"Error writing output file '{out_file}': {e}")

    print(f"✅ Saved merged file to: {out_file}")
    print(f"Join type: {args.how} | Validation: {args.validate}")
    print(f"Left keys: {left_keys} | Right keys: {right_keys}")
    print(f"Appended columns from RIGHT: {right_cols_to_add}")


if __name__ == "__main__":
    main()
```
CLI help
```txt
--left <filename>
Left/base TSV (the file to be enriched).


--right <filename>
Right/metadata TSV (the file providing columns to append).


--out <filename> (optional)
Output TSV filename. Default: <left_basename>_merged.tsv.


--left-keys <cols>
Comma‑separated key columns in the left file (e.g., Cluster or Tissue,Cluster).


--right-keys <cols>
Comma‑separated key columns in the right file. Must align one‑to‑one with --left-keys.


--map "<left:right,...>"
Alternative to --left-keys/--right-keys. Map differing key names across files
(e.g., "cluster_id:Cluster,tissue_name:Tissue").


--right-cols <cols> (optional)
Comma‑separated columns from the right file to append.
Default: all non‑key columns from the right file.


--how <left|inner|right|outer> (optional)
Join type. Default: left.


--validate <m:1|1:1|m:m> (optional)
Merge validation rule. Default: m:1.


--sep <delimiter> (optional)
Field separator for input/output. Default: tab (\t).
```
Example usage

1) Simple same‑name key:
```bash
python merge_tsv_by_keys.py \
  --left rna_single_cell_cluster.tsv \
  --right rna_single_cell_clusters.tsv \
  --left-keys Cluster \
  --right-keys Cluster \
  --right-cols "Cell count,Included in aggregation,Annotation reliability" \
  --out enriched.tsv
```

2) Multi‑column keys:
```bash
python merge_tsv_by_keys.py \
  --left fileA.tsv --right fileB.tsv \
  --left-keys "Tissue,Cluster" \
  --right-keys "Tissue,Cluster" \
  --out merged.tsv
```
3) Different key names between files:
```bash

python merge_tsv_by_keys.py \
  --left left.tsv --right right.tsv \
  --map "cluster_id:Cluster,tissue_name:Tissue" \
  --right-cols "Cell count,Annotation reliability" \
  --out merged.tsv
```

#*I Used the following command*
```bash
python merge_tsv_by_keys.py \
  --left rna_single_cell_cluster.tsv \
  --right rna_single_cell_clusters.tsv \
  --left-keys "Cluster","Tissue","Cell type" \
  --right-keys "Cluster","Tissue","Cell type" \
  --right-cols "Cell count,Included in aggregation,Annotation reliability" \
  --out  combined_expression_data.tsv
```
This gives you a file that looks like following
```tsv
Gene    Gene name       Tissue  Cluster Cell type       Read count      nCPM    Cell count      Included in aggregation Annotation reliability
ENSG00000000003 TSPAN6  ovary   c-0     ovarian stromal cells   493     92.5    8942    yes     high
ENSG00000000003 TSPAN6  ovary   c-1     ovarian stromal cells   529     80.5    6996    yes     high
ENSG00000000003 TSPAN6  ovary   c-2     ovarian stromal cells   143     52.3    6993    yes     high
ENSG00000000003 TSPAN6  ovary   c-3     ovarian stromal cells   456     91.4    6032    yes     high
ENSG00000000003 TSPAN6  ovary   c-4     vascular endothelial cells      164     28.6    5837    yes     high
```
# 2) Filtering combined dataset

Then I wanted to filter data that are not reliable. 

Human Protein atlas explains their filtration procedure as following

"Excluded from the cross-dataset aggregation and subsequent gene classification were clusters with mixed cell types, clusters with low cell type annotation confidence, and cell types within a tissue that comprised less than 30 cells or their aggregated profile contained fewer than 10,000 detected genes. We retained, however, a small number of clusters below the 30-cell threshold, provided they demonstrated more than 10,000 detected genes, to preserve representation of rare cell types. A total of 161 clusters out of 1175 clusters were excluded from the cross dataset integration and downstream analysis."

Just dropping rows with "Included in aggregation"==no (With -R "Included in aggregation" \ --keep-reliability-values "yes" \) would be enough as HPA has already done this. But I am enforcing other steps as well Just in case)

Following script tries to replicate it 

filter_integration_long.py
```py


#!/usr/bin/env python3
"""
Filter long-format expression data by cluster-level rules:
  - Drop clusters labeled as 'mixed' (substring match in Cell type; only if --drop-mixed).
  - Drop clusters with cells < min_cells OR detected genes < min_genes.
  - Exception: retain clusters below min_cells if genes >= min_genes (rare-cell preservation).
  - Reliability filter: keep or drop clusters based on values from a chosen column.
  - Included filter: row-level keep specific values; cluster-level drop clusters if any/all rows match disallowed values.

Detected genes modes:
  - compute: count unique genes per cluster with detection metric above a threshold (e.g., nCPM > 0).
  - column: use a genes count column present in the input file (per-row; aggregated per cluster).
  - map:   join a separate cluster-level mapping file that has the genes count.

QA outputs:
  - --write-summary: per-cluster table (mixed, cells, genes, reliability, included flags, keep).
  - --write-dropped: dropped clusters with explicit reasons.
"""

import argparse
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter long-format per-gene data using cluster-level rules (mixed label, min cells, min genes, reliability, 'included'), and export dropped clusters with reasons."
    )
    # I/O and format
    p.add_argument("--input", "-i", required=True, help="Path to input TSV/CSV file (long-format per gene per cluster).")
    p.add_argument("--output", "-o", required=True, help="Path to output file (filtered).")
    p.add_argument("--delimiter", "-d", default="\t", help=r"Field delimiter (default: tab '\t').")
    p.add_argument("--encoding", default="utf-8", help="File encoding (default: utf-8).")
    p.add_argument("--quotechar", default='"', help='Quote character (default: ").')
    p.add_argument("--escapechar", default=None, help="Escape character (optional).")
    p.add_argument("--verbose", "-v", action="store_true", help="Print summary details.")
    p.add_argument("--write-summary", help="Optional path to write a cluster-level summary table (CSV/TSV).")
    p.add_argument("--write-dropped", help="Optional path to write a list of dropped clusters with reasons (CSV/TSV).")

    # Columns in long-format file
    p.add_argument("--cluster-cols", nargs="+", default=["Tissue", "Cluster"],
                   help="Column(s) that identify a cluster (default: Tissue Cluster).")
    p.add_argument("--cell-type-column", default="Cell type", help="Column with cell type labels.")
    p.add_argument("--count-column", default="Cell count", help="Column with per-cluster cell count (repeated per row).")
    p.add_argument("--gene-id-column", default="Gene", help="Column with gene identifiers (used when computing detected genes).")

    # Mixed cluster filtering
    p.add_argument("--drop-mixed", action="store_true", help="Drop clusters whose cell type contains the mixed token.")
    p.add_argument("--mixed-token", default="mixed", help="Token indicating mixed cell types (default: 'mixed').")
    p.add_argument("--ignore-case", action="store_true", help="Case-insensitive matching for mixed token and reliability/included values.")

    # Thresholds & rare-cell exception
    p.add_argument("--min-cells", type=int, default=30, help="Minimum cells per cluster (default: 30).")
    p.add_argument("--min-genes", type=int, default=10000, help="Minimum detected genes per cluster (default: 10000).")
    p.add_argument("--disable-rare-exception", action="store_true",
                   help="Disable rare-cell exception (keep <min cells if genes >= min). Enabled by default.")

    # Reliability filtering (cluster-level)
    p.add_argument("--reliability-column", "-R", default="Annotation reliability",
                   help="Column to use for reliability filtering (default: 'Annotation reliability').")
    rel_group = p.add_mutually_exclusive_group()
    rel_group.add_argument("--keep-reliability-values", nargs="+",
                           help="Keep clusters whose reliability equals any of these values (e.g., 'high', 'medium high', 'medium low').")
    rel_group.add_argument("--drop-reliability-values", nargs="+",
                           help="Drop clusters whose reliability equals any of these values (e.g., 'low').")
    p.add_argument("--reliability-na-action", choices=["drop", "keep"], default="drop",
                   help="How to treat missing reliability values at cluster level (default: drop).")

    # Included filtering
    p.add_argument("--included-column", "-I", default="Included in aggregation",
                   help="Column name for 'included' status (default: 'Included in aggregation').")
    # Row-level filter (optional): keep only specified included values
    p.add_argument("--filter-included-rows", action="store_true",
                   help="Filter rows to keep only values in --included-keep-values (row-level).")
    p.add_argument("--included-keep-values", nargs="+",
                   help="Values considered 'included' when --filter-included-rows is set (e.g., yes).")
    # Back-compat convenience: treat --require-included like row-level keep==yes
    p.add_argument("--require-included", action="store_true",
                   help="(Convenience) Same as --filter-included-rows --included-keep-values yes.")
    # Cluster-level drop (optional): drop cluster if included column matches these values
    p.add_argument("--drop-included-values", nargs="+",
                   help="Drop clusters whose included column equals any of these values (e.g., 'no').")
    p.add_argument("--included-cluster-mode", choices=["any", "all"], default="any",
                   help="Drop cluster if ANY (default) or ALL rows match --drop-included-values.")

    # Detected genes modes
    p.add_argument("--genes-mode", choices=["compute", "column", "map"], default="compute",
                   help="How to obtain detected genes per cluster.")
    # compute mode
    p.add_argument("--detect-column", default="nCPM", help="Column used to decide if a gene is detected (default: 'nCPM').")
    p.add_argument("--detect-operator", choices=[">", ">="], default=">",
                   help="Operator for detection threshold (default: '>').")
    p.add_argument("--detect-threshold", type=float, default=0.0, help="Detection threshold (default: 0.0).")
    # column mode
    p.add_argument("--genes-count-column", help="Column in input that contains per-cluster detected genes (same value per cluster).")
    # map mode
    p.add_argument("--genes-mapping-file", help="Path to a separate mapping file with detected genes per cluster.")
    p.add_argument("--map-delimiter", default="\t", help="Delimiter for mapping file (default: tab).")
    p.add_argument("--map-encoding", default="utf-8", help="Encoding for mapping file (default: utf-8).")
    p.add_argument("--map-genes-count-column", help="Column in mapping file with detected genes per cluster.")

    return p.parse_args()


def _first_valid_str(series, ignore_case=False):
    """Return first non-NA string from series, normalized by stripping spaces and optional lowercasing."""
    s = series.astype("string").dropna()
    if s.size == 0:
        return pd.NA
    val = s.iloc[0]
    if pd.isna(val):
        return pd.NA
    norm = str(val).strip()
    return norm.lower() if ignore_case else norm


def _normalize_values(values, ignore_case=False):
    """Normalize a list of CLI values: strip spaces and optional lowercasing."""
    if not values:
        return []
    out = []
    for v in values:
        if v is None:
            continue
        t = str(v).strip()
        out.append(t.lower() if ignore_case else t)
    return out


def main():
    args = parse_args()

    # Read input
    try:
        df = pd.read_csv(
            args.input,
            sep=args.delimiter,
            encoding=args.encoding,
            quotechar=args.quotechar,
            escapechar=args.escapechar,
            dtype="object",
            na_filter=True
        )
    except Exception as e:
        print(f"ERROR: Failed to read input file '{args.input}': {e}", file=sys.stderr)
        sys.exit(1)

    # Translate --require-included to row-level filter (keep yes)
    if args.require_included and not args.filter_included_rows:
        args.filter_included_rows = True
        if not args.included_keep_values:
            args.included_keep_values = ["yes"]

    # Validate required columns
    missing_cols = [c for c in args.cluster_cols if c not in df.columns]
    for c in [args.cell_type_column, args.count_column]:
        if c not in df.columns:
            missing_cols.append(c)
    # Reliability column only required if a reliability filter is requested
    if (args.keep_reliability_values or args.drop_reliability_values) and (args.reliability_column not in df.columns):
        missing_cols.append(args.reliability_column)
    # Included column required if any included filter is requested
    if (args.filter_included_rows or args.drop_included_values or args.require_included) and (args.included_column not in df.columns):
        missing_cols.append(args.included_column)
    # Genes compute mode requirements
    if args.genes_mode == "compute":
        if args.gene_id_column not in df.columns:
            missing_cols.append(args.gene_id_column)
        if args.detect_column not in df.columns:
            missing_cols.append(args.detect_column)

    if missing_cols:
        print("ERROR: Missing required column(s): " + ", ".join(sorted(set(missing_cols))), file=sys.stderr)
        print("Available columns:\n  - " + "\n  - ".join(df.columns), file=sys.stderr)
        sys.exit(2)

    # --- Row-level included filter (optional) ---
    if args.filter_included_rows:
        keep_vals = _normalize_values(args.included_keep_values, ignore_case=args.ignore_case)
        included_norm = df[args.included_column].astype("string").str.strip()
        if args.ignore_case:
            included_norm = included_norm.str.lower()
        df = df[included_norm.isin(keep_vals)]

    # --- helper columns for mixed detection ---
    df["_cell_type_str"] = df[args.cell_type_column].astype("string")
    df["_is_mixed"] = df["_cell_type_str"].str.contains(
        args.mixed_token, case=not args.ignore_case, na=False, regex=False
    )

    # Group by cluster keys
    grp = df.groupby(args.cluster_cols, dropna=False)

    # Mixed flag per cluster
    mixed_flag = grp["_is_mixed"].any()
    effective_mixed_flag = mixed_flag if args.drop_mixed else mixed_flag.copy().astype(bool)
    if not args.drop_mixed:
        effective_mixed_flag[:] = False

    # Cell count per cluster: take first valid numeric value within each cluster
    cell_count_per_cluster = grp[args.count_column].agg(
        lambda s: pd.to_numeric(s, errors="coerce").dropna().iloc[0] if s.dropna().size > 0 else pd.NA
    )

    # Reliability per cluster (first non-NA, normalized)
    if args.reliability_column in df.columns:
        reliability_per_cluster = grp[args.reliability_column].agg(
            lambda s: _first_valid_str(s, ignore_case=args.ignore_case)
        )
    else:
        reliability_per_cluster = pd.Series(pd.NA, index=grp.size().index)

    # Detected genes per cluster
    if args.genes_mode == "compute":
        detect_num = pd.to_numeric(df[args.detect_column], errors="coerce")
        df["_detected"] = (detect_num >= args.detect_threshold) if args.detect_operator == ">=" else (detect_num > args.detect_threshold)
        detected_rows = df[df["_detected"].fillna(False)]
        genes_count_per_cluster = detected_rows.groupby(args.cluster_cols, dropna=False)[args.gene_id_column].nunique()
    elif args.genes_mode == "column":
        if not args.genes_count_column or args.genes_count_column not in df.columns:
            print("ERROR: --genes-count-column must be provided and exist in input for genes-mode=column.", file=sys.stderr)
            sys.exit(2)
        genes_count_per_cluster = grp[args.genes_count_column].agg(
            lambda s: pd.to_numeric(s, errors="coerce").dropna().iloc[0] if s.dropna().size > 0 else pd.NA
        )
    else:  # map mode
        if not args.genes_mapping_file or not args.map_genes_count_column:
            print("ERROR: Provide --genes-mapping-file and --map-genes-count-column for genes-mode=map.", file=sys.stderr)
            sys.exit(2)
        try:
            map_df = pd.read_csv(
                args.genes_mapping_file,
                sep=args.map_delimiter,
                encoding=args.map_encoding,
                dtype="object",
                na_filter=True
            )
        except Exception as e:
            print(f"ERROR: Failed to read mapping file '{args.genes_mapping_file}': {e}", file=sys.stderr)
            sys.exit(1)
        missing_map_keys = [c for c in args.cluster_cols if c not in map_df.columns]
        if missing_map_keys:
            print("ERROR: Mapping file missing cluster key column(s): " + ", ".join(missing_map_keys), file=sys.stderr)
            print("Mapping columns:\n  - " + "\n  - ".join(map_df.columns), file=sys.stderr)
            sys.exit(2)
        if args.map_genes_count_column not in map_df.columns:
            print(f"ERROR: Mapping file missing genes count column '{args.map_genes_count_column}'.", file=sys.stderr)
            sys.exit(2)
        map_df[args.map_genes_count_column] = pd.to_numeric(map_df[args.map_genes_count_column], errors="coerce")
        genes_count_per_cluster = map_df.set_index(args.cluster_cols)[args.map_genes_count_column]

    # Included cluster-level flags (optional) — FIXED: use helper column for boolean
    if args.drop_included_values and args.included_column in df.columns:
        included_norm = df[args.included_column].astype("string").str.strip()
        if args.ignore_case:
            included_norm = included_norm.str.lower()
        drop_vals = set(_normalize_values(args.drop_included_values, ignore_case=args.ignore_case))
        df["_included_is_drop"] = included_norm.isin(drop_vals)

        # Any / All rows in a cluster matching drop values (aggregate the helper column)
        included_bad_any = df.groupby(args.cluster_cols, dropna=False)["_included_is_drop"].any()
        included_bad_all = df.groupby(args.cluster_cols, dropna=False)["_included_is_drop"].all()
        included_bad_cluster = included_bad_any if args.included_cluster_mode == "any" else included_bad_all

        # Missing included per cluster (no non-NA value observed) — use agg to avoid FutureWarning
        has_non_na_included = grp[args.included_column].agg(lambda s: s.dropna().size > 0)
        included_missing_cluster = ~has_non_na_included
    else:
        included_bad_cluster = pd.Series(False, index=grp.size().index)
        included_missing_cluster = pd.Series(False, index=grp.size().index)

    # Align all series to a common MultiIndex covering all clusters in the input
    cluster_index = grp.size().index
    summary = pd.DataFrame({
        "mixed_flag": mixed_flag.reindex(cluster_index, fill_value=False).astype(bool),
        "effective_mixed_flag": effective_mixed_flag.reindex(cluster_index, fill_value=False).astype(bool),
        "cell_count": pd.to_numeric(cell_count_per_cluster.reindex(cluster_index), errors="coerce"),
        "genes_count": pd.to_numeric(genes_count_per_cluster.reindex(cluster_index), errors="coerce"),
        "reliability": reliability_per_cluster.reindex(cluster_index),
        "included_bad": included_bad_cluster.reindex(cluster_index, fill_value=False).astype(bool),
        "included_missing": included_missing_cluster.reindex(cluster_index, fill_value=False).astype(bool),
    })
    summary.index.names = args.cluster_cols

    # Apply thresholds (treat NA as failing/too few unless exception applies)
    too_few_cells = (summary["cell_count"] < args.min_cells)
    too_few_genes = (summary["genes_count"] < args.min_genes)
    too_few_cells_fill = too_few_cells.fillna(True)
    too_few_genes_fill = too_few_genes.fillna(True)

    # Rare-cell exception: keep clusters with <min cells but >=min genes
    rare_exception_on = not args.disable_rare_exception
    rare_keep = (too_few_cells_fill) & (~too_few_genes_fill) if rare_exception_on else pd.Series(False, index=summary.index)

    # Reliability filter
    reliability_ok = pd.Series(True, index=summary.index)
    if args.keep_reliability_values or args.drop_reliability_values:
        allowed = _normalize_values(args.keep_reliability_values, ignore_case=args.ignore_case) if args.keep_reliability_values else None
        disallowed = _normalize_values(args.drop_reliability_values, ignore_case=args.ignore_case) if args.drop_reliability_values else None
        rel_norm = summary["reliability"].astype("string").str.strip()
        if args.ignore_case:
            rel_norm = rel_norm.str.lower()
        if allowed is not None:
            reliability_ok = rel_norm.isin(allowed)
        else:
            reliability_ok = ~rel_norm.isin(disallowed)
        reliability_ok = reliability_ok.fillna(True if args.reliability_na_action == "keep" else False)

    # Final keep decision (include cluster-level 'included' drop)
    keep_clusters = reliability_ok & (~summary["effective_mixed_flag"]) & (~summary["included_bad"]) & (
        ((~too_few_cells_fill) & (~too_few_genes_fill)) | rare_keep
    )
    summary["keep"] = keep_clusters

    # For QA: explicit reasons for dropped clusters
    dropped = summary[~summary["keep"]].copy()
    if not dropped.empty:
        dropped["missing_cells"] = summary["cell_count"].isna().reindex(dropped.index)
        dropped["missing_genes"] = summary["genes_count"].isna().reindex(dropped.index)
        dropped["low_cells"] = too_few_cells.reindex(dropped.index).fillna(False)
        dropped["low_genes"] = too_few_genes.reindex(dropped.index).fillna(False)
        dropped["mixed"] = summary["effective_mixed_flag"].reindex(dropped.index).fillna(False)
        dropped["missing_reliability"] = dropped["reliability"].isna()
        dropped["bad_included"] = summary["included_bad"].reindex(dropped.index).fillna(False)
        dropped["missing_included"] = summary["included_missing"].reindex(dropped.index).fillna(False)

        # Recompute reliability_ok on dropped index to get a 'bad_reliability' flag
        rel_norm_dropped = dropped["reliability"].astype("string").str.strip()
        if args.ignore_case:
            rel_norm_dropped = rel_norm_dropped.str.lower()
        if args.keep_reliability_values:
            allowed = _normalize_values(args.keep_reliability_values, ignore_case=args.ignore_case)
            bad_rel = ~rel_norm_dropped.isin(allowed)
        elif args.drop_reliability_values:
            disallowed = _normalize_values(args.drop_reliability_values, ignore_case=args.ignore_case)
            bad_rel = rel_norm_dropped.isin(disallowed)
        else:
            bad_rel = pd.Series(False, index=dropped.index)
        bad_rel = bad_rel.fillna(args.reliability_na_action == "drop")
        dropped["bad_reliability"] = bad_rel

        def _combine_reasons(row):
            reasons = []
            if row.get("mixed", False): reasons.append("mixed")
            if row.get("bad_included", False): reasons.append("bad_included")
            elif row.get("missing_included", False): reasons.append("missing_included")
            if row.get("bad_reliability", False): reasons.append("bad_reliability")
            elif row.get("missing_reliability", False): reasons.append("missing_reliability")
            if row.get("low_cells", False): reasons.append("low_cells")
            elif row.get("missing_cells", False): reasons.append("missing_cells")
            if row.get("low_genes", False): reasons.append("low_genes")
            elif row.get("missing_genes", False): reasons.append("missing_genes")
            return "; ".join(reasons) if reasons else "unspecified"

        dropped["reason"] = dropped.apply(_combine_reasons, axis=1)

    # Filter long-format rows by cluster membership
    df_idx = df.set_index(args.cluster_cols)
    kept_idx = summary.index[summary["keep"]]
    filtered = df_idx.loc[df_idx.index.isin(kept_idx)].reset_index()

    # Drop helper columns before writing
    for col in ["_cell_type_str", "_is_mixed", "_detected", "_included_is_drop"]:
        if col in filtered.columns:
            filtered = filtered.drop(columns=[col])

    if args.verbose:
        total_clusters = len(summary)
        kept_clusters = int(summary["keep"].sum())
        dropped_clusters = total_clusters - kept_clusters
        total_rows = len(df)
        kept_rows = len(filtered)
        dropped_rows = total_rows - kept_rows
        rel_mode = ("KEEP " + ", ".join(args.keep_reliability_values)) if args.keep_reliability_values \
                   else (("DROP " + ", ".join(args.drop_reliability_values)) if args.drop_reliability_values else "OFF")
        inc_mode = (f"DROP values={args.drop_included_values} mode={args.included_cluster_mode}") if args.drop_included_values else "OFF"
        row_inc = (f"ROW keep={args.included_keep_values}") if args.filter_included_rows else "OFF"
        print(f"[INFO] Clusters total: {total_clusters}, kept: {kept_clusters}, dropped: {dropped_clusters}")
        print(f"[INFO] Mixed filtering: {'ON' if args.drop_mixed else 'OFF'} (token='{args.mixed_token}', ignore_case={args.ignore_case})")
        print(f"[INFO] Thresholds: min_cells={args.min_cells}, min_genes={args.min_genes}")
        print(f"[INFO] Reliability column: {args.reliability_column}; filter: {rel_mode}; NA policy={args.reliability_na_action}")
        print(f"[INFO] Included column: {args.included_column}; cluster filter: {inc_mode}; row filter: {row_inc}")
        print(f"[INFO] Rare-cell exception: {'ON' if not args.disable_rare_exception else 'OFF'}")
        print(f"[INFO] Rows total: {total_rows}, kept: {kept_rows}, dropped: {dropped_rows}")

    # Optional summary output
    if args.write_summary:
        try:
            sep = args.delimiter
            summary.reset_index().to_csv(args.write_summary, sep=sep, encoding=args.encoding, index=False)
            if args.verbose:
                print(f"[INFO] Wrote cluster summary to '{args.write_summary}'")
        except Exception as e:
            print(f"ERROR: Failed to write cluster summary '{args.write_summary}': {e}", file=sys.stderr)

    # Optional dropped clusters output
    if args.write_dropped:
        try:
            sep = args.delimiter
            out = dropped.reset_index()[args.cluster_cols + [
                "mixed_flag", "effective_mixed_flag", "cell_count", "genes_count",
                "reliability", "included_bad", "included_missing",
                "mixed", "low_cells", "low_genes",
                "missing_cells", "missing_genes", "missing_reliability",
                "bad_reliability", "bad_included", "missing_included", "reason"
            ]]
        except KeyError:
            # If 'dropped' is empty, still create an empty table with the expected columns
            out_cols = args.cluster_cols + [
                "mixed_flag", "effective_mixed_flag", "cell_count", "genes_count",
                "reliability", "included_bad", "included_missing",
                "mixed", "low_cells", "low_genes",
                "missing_cells", "missing_genes", "missing_reliability",
                "bad_reliability", "bad_included", "missing_included", "reason"
            ]
            out = pd.DataFrame(columns=out_cols)
        try:
            out.to_csv(args.write_dropped, sep=sep, encoding=args.encoding, index=False)
            if args.verbose:
                print(f"[INFO] Wrote dropped clusters list to '{args.write_dropped}'")
        except Exception as e:
            print(f"ERROR: Failed to write dropped clusters '{args.write_dropped}': {e}", file=sys.stderr)

    # Write filtered long-format data
    try:
        filtered.to_csv(
            args.output,
            sep=args.delimiter,
            encoding=args.encoding,
            index=False,
            quoting=0  # csv.QUOTE_MINIMAL
        )
        if args.verbose:
            print(f"[INFO] Wrote {len(filtered)} rows to '{args.output}'")
    except Exception as e:
        print(f"ERROR: Failed to write output file '{args.output}': {e}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()

```
CLI help
```txt

usage: filter_integration_long.py [-h] --input INPUT --output OUTPUT
                                  [--delimiter DELIMITER] [--encoding ENCODING] [--quotechar QUOTECHAR]
                                  [--escapechar ESCAPECHAR] [--verbose]
                                  [--write-summary WRITE_SUMMARY] [--write-dropped WRITE_DROPPED]
                                  [--cluster-cols CLUSTER_COLS [CLUSTER_COLS ...]]
                                  [--cell-type-column CELL_TYPE_COLUMN] [--count-column COUNT_COLUMN]
                                  [--gene-id-column GENE_ID_COLUMN] [--drop-mixed] [--mixed-token MIXED_TOKEN]
                                  [--ignore-case] [--min-cells MIN_CELLS] [--min-genes MIN_GENES]
                                  [--disable-rare-exception] [--reliability-column RELIABILITY_COLUMN]
                                  [--keep-reliability-values KEEP_RELIABILITY_VALUES [KEEP_RELIABILITY_VALUES ...] |
                                   --drop-reliability-values DROP_RELIABILITY_VALUES [DROP_RELIABILITY_VALUES ...]]
                                  [--reliability-na-action {drop,keep}]
                                  [--included-column INCLUDED_COLUMN] [--filter-included-rows]
                                  [--included-keep-values INCLUDED_KEEP_VALUES [INCLUDED_KEEP_VALUES ...]]
                                  [--require-included] [--drop-included-values DROP_INCLUDED_VALUES [DROP_INCLUDED_VALUES ...]]
                                  [--included-cluster-mode {any,all}]
                                  [--genes-mode {compute,column,map}] [--detect-column DETECT_COLUMN]
                                  [--detect-operator {>,>=}] [--detect-threshold DETECT_THRESHOLD]
                                  [--genes-count-column GENES_COUNT_COLUMN]
                                  [--genes-mapping-file GENES_MAPPING_FILE] [--map-delimiter MAP_DELIMITER]
                                  [--map-encoding MAP_ENCODING] [--map-genes-count-column MAP_GENES_COUNT_COLUMN]

Filter long-format per-gene data using cluster-level rules (mixed label, min cells, min genes,
reliability, and 'included'), with QA outputs.

required arguments:
  --input, -i                Path to input TSV/CSV (long-format per gene per cluster).
  --output, -o               Path to output TSV/CSV (filtered).

format & IO:
  --delimiter, -d            Field delimiter (default: tab '\t').
  --encoding                 File encoding (default: utf-8).
  --quotechar                Quote character (default: ").
  --escapechar               Escape character (optional).
  --verbose, -v              Print summary details to stdout.
  --write-summary            Path to save a per-cluster summary table.
  --write-dropped            Path to save dropped clusters with explicit reasons.

columns (input schema):
  --cluster-cols             Column(s) identifying a cluster (default: Tissue Cluster).
                             Pass multiple columns separated by spaces.
  --cell-type-column         Column with cell type labels (default: "Cell type").
  --count-column             Column with per-cluster cell count (default: "Cell count").
  --gene-id-column           Column with gene identifiers (default: "Gene").

mixed cluster filtering:
  --drop-mixed               Drop clusters whose cell type contains the mixed token.
  --mixed-token              Token indicating mixed clusters (default: "mixed").
  --ignore-case              Case-insensitive matching for mixed token, reliability, and included values.

thresholds & rare-cell exception:
  --min-cells                Minimum cells per cluster (default: 30).
  --min-genes                Minimum detected genes per cluster (default: 10000).
  --disable-rare-exception   Turn OFF the rare-cell exception (ON by default).
                             Rare-cell exception: keep clusters with <min cells if genes ≥ min.

reliability filtering (cluster-level):
  --reliability-column, -R   Column used for reliability (default: "Annotation reliability").
  --keep-reliability-values  Keep only clusters whose reliability equals any of these values
                             (e.g., high "medium high" "medium low"). (mutually exclusive with --drop-reliability-values)
  --drop-reliability-values  Drop clusters whose reliability equals any of these values
                             (e.g., low). (mutually exclusive with --keep-reliability-values)
  --reliability-na-action    How to treat missing reliability at cluster level: drop | keep (default: drop).

included filtering (row- and cluster-level):
  --included-column, -I      Column used for included status (default: "Included in aggregation").
  --filter-included-rows     Row-level filter: keep only rows whose included value is in --included-keep-values.
  --included-keep-values     Values considered included for row-level filtering (e.g., yes).
  --require-included         Convenience: same as "--filter-included-rows --included-keep-values yes".
  --drop-included-values     Cluster-level filter: drop clusters if the included value matches any of these (e.g., no).
  --included-cluster-mode    When using --drop-included-values, drop the cluster if ANY row matches (default: any),
                             or only if ALL rows match (all).

detected genes (choose one mode):
  --genes-mode               How to obtain detected genes per cluster: compute | column | map (default: compute).

  compute mode:
    --detect-column          Column used for detection metric (default: nCPM).
    --detect-operator        Detection operator: > | >= (default: >).
    --detect-threshold       Detection threshold (float, default: 0.0). Example: nCPM >= 1.

  column mode:
    --genes-count-column     Column in input containing per-cluster detected gene counts.

  map mode:
    --genes-mapping-file     Separate mapping file with detected genes per cluster.
    --map-delimiter          Mapping file delimiter (default: tab).
    --map-encoding           Mapping file encoding (default: utf-8).
    --map-genes-count-column Column in mapping file with detected gene counts.

notes:
  • Cluster identity = unique combination of the columns passed to --cluster-cols (e.g., Tissue + Cluster).
  • Values with spaces MUST be quoted in the shell: "Included in aggregation", "medium high".
  • Missing values in cell count / detected genes are treated as failing thresholds (i.e., dropped),
    unless preserved by the rare-cell exception (cells < min AND genes ≥ min).
  • Row-level included filtering (--filter-included-rows) removes rows before summarizing clusters.
    Cluster-level included filtering (--drop-included-values) records QA reasons (bad_included / missing_included).

examples:
  # Typical use: compute detected genes (nCPM > 0), drop mixed, drop clusters with any 'Included in aggregation' == no
  python filter_integration_long.py \
    --input combined_expression_data.tsv \
    --output integrated_filtered.tsv \
    --cluster-cols Tissue Cluster \
    --min-cells 30 --min-genes 10000 \
    --drop-mixed --ignore-case \
    --included-column "Included in aggregation" \
    --drop-included-values no \
    --included-cluster-mode any \
    --write-summary cluster_summary.tsv \
    --write-dropped dropped_clusters.tsv \
    --verbose

  # Enforce reliability: keep only high + medium high + medium low; compute detected genes with nCPM >= 1
  python filter_integration_long.py \
    --input combined_expression_data.tsv \
    --output integrated_filtered.tsv \
    --cluster-cols Tissue Cluster \
    --min-cells 30 --min-genes 10000 \
    --drop-mixed --ignore-case \
    -R "Annotation reliability" \
    --keep-reliability-values high "medium high" "medium low" \
    --genes-mode compute --detect-column nCPM --detect-operator '>=' --detect-threshold 1 \
    --write-summary cluster_summary.tsv \
    --write-dropped dropped_clusters.tsv \
    --verbose

  # Use a precomputed per-cluster genes count in the input (column mode)
  python filter_integration_long.py \
    --input combined_expression_data.tsv \
    --output integrated_filtered.tsv \
    --cluster-cols Tissue Cluster \
    --genes-mode column --genes-count-column "Detected genes" \
    --min-cells 30 --min-genes 10000 \
    --drop-mixed \
    --verbose

```
#*I am running it like following*
```bash
python filter_integration_long.py \
  --input combined_expression_data.tsv \
  --output integrated_filtered.tsv \
  --cluster-cols Tissue Cluster \
  --min-cells 30 --min-genes 10000 \
  --drop-mixed --ignore-case \
  --included-column "Annotation reliability" \
  --drop-included-values "low" "medium low" \
  --included-cluster-mode any \
  --write-summary cluster_summary.tsv \
  --write-dropped dropped_clusters.tsv \
  --verbose
```
This dropped a total of 161 clusters out of 1175 clusters exactly like HPA pipeline had done

# 3) Calculations

1) First, with the resulting dataset from the previous step containing multiple clusters for each gene-cell type combination, I wanted to caclulate weighted mean by cell count number (Just like what HPA did) for each Tissue*Cell type*Gene combination

I did that with aggregate_within_dataset.py

aggregate_within_dataset.py
```py

#!/usr/bin/env python3
"""
Within-Dataset aggregation (weighted mean by cell count):

For each (dataset, cell type, gene), compute for one or more value columns:
    aggregated_<col> = sum_j( Cell count_j * <col>_j ) / sum_j( Cell count_j )

Defaults assume:
  - Dataset: "Tissue"
  - Cell type: "Cell type"
  - Gene ID: "Gene"
  - Values: ["nCPM"] (change with --agg-values, e.g., include "Read count")
  - Weight: "Cell count"

Features:
  - Filter rows by "Included in aggregation" values (e.g., keep only 'yes').
  - Case-insensitive comparisons (for included filter).
  - NA handling policies for values and weights.
  - Optional stats (sum of weights, number of clusters).
  - Optional pivot to a wide matrix (genes/cell types × dataset × metric).
  - Include extra passthrough columns (e.g., "Gene name") via --extra-cols.
"""

import argparse
import sys
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute weighted mean per (dataset, cell type, gene) using Cell count as weights for one or more value columns."
    )
    # I/O
    p.add_argument("--input", "-i", required=True, help="Path to input TSV/CSV (long-format per gene per cluster).")
    p.add_argument("--output", "-o", required=True, help="Path to output TSV/CSV (aggregated).")
    p.add_argument("--delimiter", "-d", default="\t", help=r"Field delimiter (default: tab '\t').")
    p.add_argument("--encoding", default="utf-8", help="File encoding (default: utf-8).")
    p.add_argument("--quotechar", default='"', help='Quote character (default: ").')
    p.add_argument("--escapechar", default=None, help="Escape character (optional).")
    p.add_argument("--verbose", "-v", action="store_true", help="Print summary details.")

    # Columns
    p.add_argument("--dataset-col", default="Tissue", help="Column representing dataset (default: 'Tissue').")
    p.add_argument("--cell-type-col", default="Cell type", help="Column for cell type (default: 'Cell type').")
    p.add_argument("--gene-id-col", default="Gene", help="Column for gene IDs (default: 'Gene').")
    p.add_argument("--weight-col", default="Cell count", help="Column for weights (default: 'Cell count').")
    p.add_argument("--cluster-col", default="Cluster", help="Optional column for cluster ID (used in stats).")

    # Values to aggregate (one or more)
    p.add_argument("--agg-values", nargs="+", default=["nCPM"],
                   help="One or more value columns to aggregate by weighted mean (e.g., nCPM \"Read count\").")

    # NEW: extra passthrough columns (e.g., Gene name)
    p.add_argument("--extra-cols", nargs="+",
                   help="Additional columns to include in the output by taking the first non-NA value per group (e.g., \"Gene name\").")

    # Optional: included filter
    p.add_argument("--included-column", default="Included in aggregation",
                   help="Column name to check row inclusion (default: 'Included in aggregation').")
    p.add_argument("--require-included", nargs="+",
                   help="Keep only rows where included-column equals one of these values (e.g., yes).")
    p.add_argument("--ignore-case", action="store_true",
                   help="Case-insensitive comparisons for included values.")

    # Policies & output
    p.add_argument("--value-na-policy", choices=["drop", "zero"], default="drop",
                   help="How to handle NA in value columns during aggregation (default: drop).")
    p.add_argument("--weight-na-policy", choices=["drop", "zero"], default="drop",
                   help="How to handle NA in weight column during aggregation (default: drop).")
    p.add_argument("--add-stats", action="store_true",
                   help="Add stats columns: weight_sum and n_clusters.")
    p.add_argument("--precision", type=int, default=6, help="Decimal rounding for aggregated values (default: 6).")

    # Pivot options
    p.add_argument("--pivot", choices=["none", "wide"], default="none",
                   help="Output shape: none (long) or wide (pivot).")
    p.add_argument("--pivot-axis", choices=["genes", "cell_types"], default="genes",
                   help="For wide pivot: rows are genes (default) or cell types.")

    return p.parse_args()


def _normalize_series_text(s: pd.Series, lower: bool):
    s = s.astype("string").str.strip()
    return s.str.lower() if lower else s


def main():
    args = parse_args()

    # Read input
    try:
        df = pd.read_csv(
            args.input,
            sep=args.delimiter,
            encoding=args.encoding,
            quotechar=args.quotechar,
            escapechar=args.escapechar,
            dtype="object",  # read as strings first
            na_filter=True
        )
    except Exception as e:
        print(f"ERROR: Failed to read input file '{args.input}': {e}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required = [args.dataset_col, args.cell_type_col, args.gene_id_col, args.weight_col]
    missing = [c for c in required if c not in df.columns]
    # Also require all agg-values columns
    missing += [c for c in args.agg_values if c not in df.columns]
    if args.extra_cols:
        missing += [c for c in args.extra_cols if c not in df.columns]
    if missing:
        print("ERROR: Missing required column(s): " + ", ".join(sorted(set(missing))), file=sys.stderr)
        print("Available columns:\n  - " + "\n  - ".join(df.columns), file=sys.stderr)
        sys.exit(2)

    # Optional: included filter
    if args.require_included:
        if args.included_column not in df.columns:
            print(f"ERROR: Included column '{args.included_column}' not found.", file=sys.stderr)
            sys.exit(2)
        keep_vals = [str(v).strip() for v in args.require_included]
        if args.ignore_case:
            keep_vals = [v.lower() for v in keep_vals]
        inc_norm = _normalize_series_text(df[args.included_column], lower=args.ignore_case)
        df = df[inc_norm.isin(keep_vals)]

    # Numeric casts (robust) — weights
    wt = pd.to_numeric(df[args.weight_col], errors="coerce")

    # NA policies & mask: start from weights
    mask_valid = pd.Series(True, index=df.index)
    if args.weight_na_policy == "drop":
        mask_valid &= wt.notna()
    else:  # zero
        wt = wt.fillna(0.0)

    # Numeric casts (robust) — for each value column
    val_num = {}
    for col in args.agg_values:
        vc = pd.to_numeric(df[col], errors="coerce")
        if args.value_na_policy == "drop":
            mask_valid &= vc.notna()
        else:  # zero
            vc = vc.fillna(0.0)
        val_num[col] = vc

    # Filtered DataFrame that will be used for grouping
    df_valid = df[mask_valid].copy()
    if df_valid.empty:
        print("WARN: No rows remaining after NA policies/filtering.", file=sys.stderr)

    # Attach numeric helper columns so groupby uses numeric data
    df_valid["_weight_num"] = wt.loc[df_valid.index]
    for col in args.agg_values:
        df_valid[f"_{col}_num"] = val_num[col].loc[df_valid.index]
        df_valid[f"_w_times_{col}"] = df_valid["_weight_num"] * df_valid[f"_{col}_num"]

    # Group keys
    keys = [args.dataset_col, args.cell_type_col, args.gene_id_col]

    # Build aggregation dict dynamically
    agg_dict = {"weight_sum": ("_weight_num", "sum")}
    if args.cluster_col in df_valid.columns:
        agg_dict["n_clusters"] = (args.cluster_col, "nunique")
    else:
        agg_dict["n_clusters"] = ("_weight_num", "size")  # fallback to row count

    for col in args.agg_values:
        agg_dict[f"weighted_sum_{col}"] = (f"_w_times_{col}", "sum")

    # Aggregate passthrough extra columns (first non-NA per group)
    def _first_valid(s):
        s = s.dropna()
        return s.iloc[0] if len(s) else pd.NA

    extra_results = {}
    if args.extra_cols:
        for ex in args.extra_cols:
            extra_results[ex] = (ex, _first_valid)

    grouped = df_valid.groupby(keys, dropna=False).agg(**agg_dict, **extra_results)

    # Compute weighted means safely (denom > 0 → num/denom, else NaN), and round
    aggregated = grouped.copy()
    denom = aggregated["weight_sum"]
    for col in args.agg_values:
        num = aggregated[f"weighted_sum_{col}"]
        outcol = f"aggregated_{col}"
        aggregated[outcol] = np.where(denom > 0, num / denom, np.nan)
        aggregated[outcol] = aggregated[outcol].round(args.precision)

    # Prepare LONG output by default (include extras if requested)
    out_cols = keys.copy()
    if args.extra_cols:
        out_cols += args.extra_cols  # e.g., "Gene name"
    out_cols += [f"aggregated_{col}" for col in args.agg_values]
    if args.add_stats:
        out_cols += ["weight_sum", "n_clusters"]

    result_long = aggregated.reset_index()[out_cols]

    # Optional pivot to wide format
    if args.pivot == "wide":
        # Melt metrics to a single column, then pivot
        metric_cols = [f"aggregated_{col}" for col in args.agg_values]
        id_vars = keys + (args.extra_cols if args.extra_cols else [])
        if args.add_stats:
            id_vars += ["weight_sum", "n_clusters"]
        long_melt = result_long.melt(id_vars=id_vars,
                                     value_vars=metric_cols,
                                     var_name="metric",
                                     value_name="value")
        # Simplify metric names: 'aggregated_nCPM' -> 'nCPM', etc.
        long_melt["metric"] = long_melt["metric"].str.replace(r"^aggregated_", "", regex=True)

        if args.pivot_axis == "genes":
            pivot = long_melt.pivot_table(
                index=[args.gene_id_col] + ([c for c in args.extra_cols if c != args.gene_id_col] if args.extra_cols else []),
                columns=[args.dataset_col, args.cell_type_col, "metric"],
                values="value",
                aggfunc="first"
            )
            pivot.columns = ["{}|{}|{}".format(str(ds), str(ct), str(m)) for (ds, ct, m) in pivot.columns]
            pivot = pivot.reset_index()
        else:
            pivot = long_melt.pivot_table(
                index=[args.cell_type_col],
                columns=[args.dataset_col, args.gene_id_col, "metric"],
                values="value",
                aggfunc="first"
            )
            pivot.columns = ["{}|{}|{}".format(str(ds), str(g), str(m)) for (ds, g, m) in pivot.columns]
            pivot = pivot.reset_index()

        try:
            pivot.to_csv(args.output, sep=args.delimiter, encoding=args.encoding, index=False)
            if args.verbose:
                print(f"[INFO] Wrote wide matrix to '{args.output}' with shape {pivot.shape}")
        except Exception as e:
            print(f"ERROR: Failed to write output file '{args.output}': {e}", file=sys.stderr)
            sys.exit(3)
        return

    # Write LONG table
    try:
        result_long.to_csv(args.output, sep=args.delimiter, encoding=args.encoding, index=False)
        if args.verbose:
            groups = len(aggregated)
            zero_denom = int((aggregated["weight_sum"] <= 0).sum())
            print(f"[INFO] Aggregated {groups} (dataset, cell type, gene) groups")
            if zero_denom:
                print(f"[WARN] {zero_denom} groups had weight_sum <= 0; aggregated values set to NaN.")
            print(f"[INFO] Wrote aggregated table to '{args.output}' with {len(result_long)} rows")
    except Exception as e:
        print(f"ERROR: Failed to write output file '{args.output}': {e}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
```
CLI help
```txt

usage: aggregate_within_dataset.py [-h] --input INPUT --output OUTPUT
                                   [--delimiter DELIMITER] [--encoding ENCODING]
                                   [--quotechar QUOTECHAR] [--escapechar ESCAPECHAR]
                                   [--verbose] [--dataset-col DATASET_COL]
                                   [--cell-type-col CELL_TYPE_COL] [--gene-id-col GENE_ID_COL]
                                   [--weight-col WEIGHT_COL] [--cluster-col CLUSTER_COL]
                                   [--agg-values AGG_VALUES [AGG_VALUES ...]]
                                   [--extra-cols EXTRA_COLS [EXTRA_COLS ...]]
                                   [--included-column INCLUDED_COLUMN]
                                   [--require-included REQUIRE_INCLUDED [REQUIRE_INCLUDED ...]]
                                   [--ignore-case]
                                   [--value-na-policy {drop,zero}]
                                   [--weight-na-policy {drop,zero}]
                                   [--add-stats] [--precision PRECISION]
                                   [--pivot {none,wide}] [--pivot-axis {genes,cell_types}]

Compute weighted mean per (dataset, cell type, gene) using Cell count as weights.
```
```txt
Required

--input, -i — Path to input TSV/CSV (long-format per gene per cluster).
--output, -o — Path to output TSV/CSV (aggregated).

Columns

--dataset-col — Dataset column (default: Tissue).
--cell-type-col — Cell type column (default: "Cell type").
--gene-id-col — Gene ID column (default: Gene).
--weight-col — Weight column (default: "Cell count").
--cluster-col — Cluster ID column (optional; used to count contributing clusters).

Values to aggregate

--agg-values — One or more value columns to aggregate by weighted mean (default: nCPM).
Example: --agg-values nCPM "Read count"

Extra passthrough columns

--extra-cols — Additional columns to include in the output by taking the first non‑NA per group.
Example: --extra-cols "Gene name" "Annotation reliability"

Included filter

--included-column — Column to check inclusion (default: "Included in aggregation").
--require-included — Keep only rows with these values (e.g., yes).
--ignore-case — Case-insensitive matching for included values.

NA handling & output

--value-na-policy — How to handle NA in value columns: drop (default) | zero.
--weight-na-policy — How to handle NA in weights: drop (default) | zero.
--add-stats — Add weight_sum (∑ weights) and n_clusters (distinct clusters) to output.
--precision — Decimal rounding for aggregated values (default: 6).

Format & pivot

--delimiter, -d — Field delimiter (default: tab).
--encoding — File encoding (default: UTF‑8).
--quotechar, --escapechar — CSV parsing options.
--pivot — Output shape: none (long, default) | wide (pivoted matrix).
--pivot-axis — For wide pivot, rows are genes (default) or cell_types.
--verbose, -v — Print summary details.
```
Examples

HPA‑style weighted mean of nCPM by Cell count (long output) + include Gene name:
```bash

python aggregate_within_dataset.py \
  --input integrated_filtered.tsv \
  --output aggregated_within_dataset.tsv \
  --dataset-col Tissue \
  --cell-type-col "Cell type" \
  --gene-id-col Gene \
  --weight-col "Cell count" \
  --agg-values nCPM \
  --extra-cols "Gene name" \
  --require-included yes \
  --ignore-case \
  --add-stats \
  --precision 6 \
  --verbose
```
Aggregate nCPM and Read count in one run:
```bash

python aggregate_within_dataset.py \
  --input integrated_filtered.tsv \
  --output aggregated_within_dataset.tsv \
  --dataset-col Tissue \
  --cell-type-col "Cell type" \
  --gene-id-col Gene \
  --weight-col "Cell count" \
  --agg-values nCPM "Read count" \
  --extra-cols "Gene name" \
  --require-included yes --ignore-case \
  --add-stats --precision 6 --verbose
```
Wide matrix (rows=genes, columns=dataset|cell type|metric):
```bash

python aggregate_within_dataset.py \
  --input integrated_filtered.tsv \
  --output aggregated_matrix.tsv \
  --dataset-col Tissue \
  --cell-type-col "Cell type" \
  --gene-id-col Gene \
  --weight-col "Cell count" \
  --agg-values nCPM "Read count" \
  --extra-cols "Gene name" \
  --require-included yes --ignore-case \
  --pivot wide --pivot-axis genes \
  --precision 6 --verbose
```
Notes

Weighting: All metrics in --agg-values are weighted by --weight-col (typically "Cell count").
Grouping keys: (dataset-col, cell-type-col, gene-id-col).
Zero weights: Groups with weight_sum <= 0 yield NaN for aggregated values (logged with --verbose).
Quoting: Quote argument values that contain spaces: "Cell type", "Gene name", "Included in aggregation".

#*I ran it like following*
```bash
python aggregate_within_dataset.py \
  --input integrated_filtered.tsv \
  --output aggregated_within_dataset.tsv \
  --dataset-col Tissue \
  --cell-type-col "Cell type" \
  --gene-id-col Gene \
  --weight-col "Cell count" \
  --agg-values nCPM \
  --extra-cols "Gene name" \
  --require-included yes \
  --ignore-case \
  --add-stats \
  --precision 6 \
  --verbose
```
