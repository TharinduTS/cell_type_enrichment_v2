# Cell_Type_enrichment_V_2.1 

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

## 1-I Introduction
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

## 1-II Merge script
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
## 1-III CLI help
CLI help

#**************	NOTE ************************************

This script can leave rows without matching keys on right table column value empty. Therefore you should use a command like following to find the rows with unassigned values to avoid errors down the pipeline (I have dealt with this issue in chapter 5)
see the example below from section 5

```bash
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Cell type group") col=i; next} $col=="" {print}' enrich_values_with_cell_class_data.tsv | cut -f 3 | sort -u
```
#********************************************************
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
## 1-IV Used command
#*I Used the following command*
```bash
python merge_tsv_by_keys.py \
  --left rna_single_cell_cluster.tsv \
  --right rna_single_cell_clusters.tsv \
  --left-keys "Cluster","Tissue","Cell type" \
  --right-keys "Cluster","Tissue","Cell type" \
  --right-cols "Cell type class,Cell type detail,Cell count,Annotation reliability" \
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

## 2-I Introduction
Then I wanted to filter data that are not reliable. 

Human Protein atlas explains their filtration procedure as following

"Excluded from the cross-dataset aggregation and subsequent gene classification were clusters with mixed cell types, clusters with low cell type annotation confidence, and cell types within a tissue that comprised less than 30 cells or their aggregated profile contained fewer than 10,000 detected genes. We retained, however, a small number of clusters below the 30-cell threshold, provided they demonstrated more than 10,000 detected genes, to preserve representation of rare cell types. A total of 161 clusters out of 1175 clusters were excluded from the cross dataset integration and downstream analysis."

I followed the same filtering parameters 

Following script tries to replicate thier filtering methods

## 2-II Filtering Script
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
## 2-III CLI help
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
## 2-IV Run command
#*I am running it like following*
```bash
python filter_integration_long.py \
  --input combined_expression_data.tsv \
  --output combined_expression_data_filtered.tsv \
  --cluster-cols Tissue Cluster \
  --min-cells 30 --min-genes 8000 \
  --drop-mixed --ignore-case \
  --included-column "Annotation reliability" \
  --drop-included-values "low" "medium low" \
  --included-cluster-mode any \
  --write-summary cluster_summary.tsv \
  --write-dropped dropped_clusters.tsv \
  --verbose
```
This dropped a total of 161 clusters out of 1175 clusters exactly like HPA pipeline had done when I used the same parameters (But in this version I set gene limit to 8000 to include platelets and  megakaryocytes)

# 3) Cell type enrichment

## 3-I Introduction
Then I used my script celltype_enrichment_v1_4.py to calculate weighted nCPM by cell count, celltype enrichment, select top gene-cell type combinations, enforce median scaling and safeguards

## 3-II Celltype enrichment script
celltype_enrichment_v1_4.py

```py

#!/usr/bin/env python3
"""
Cell-type–aware aggregation and enrichment with optional weighting
and Yanai's τ (specificity) — v1.4

Includes:
- τ report column (specificity_tau)
- Modes: --specificity-mode off|filter|penalize
- Safe log2 with epsilon + zero masking
- Series-based τ computation to avoid FutureWarning
"""

import argparse
import numpy as np
import pandas as pd
from typing import Tuple, List, Optional

# ---------- Utilities ----------

def coerce_numeric(series: pd.Series) -> Tuple[pd.Series, List]:
    s = pd.to_numeric(series, errors="coerce")
    non_numeric = series.loc[s.isna()].unique().tolist()
    return s, non_numeric


def drop_genes_with_no_expression(
    agg_df: pd.DataFrame,
    expr_col: Optional[str] = None,
    treat_nan_as_zero: bool = False,
) -> Tuple[pd.DataFrame, List[str], int]:
    if expr_col is None:
        expr_col = "avg_nCPM" if "avg_nCPM" in agg_df.columns else "nCPM"
    df = agg_df.copy()
    df[expr_col] = pd.to_numeric(df[expr_col], errors="coerce")
    df[expr_col] = df[expr_col].replace([np.inf, -np.inf], np.nan)
    expr_for_test = df[expr_col].fillna(0) if treat_nan_as_zero else df[expr_col]
    gene_max = expr_for_test.groupby(df["Gene"]).max()
    genes_all_zero = gene_max[gene_max == 0].index.tolist()
    before = len(df)
    filtered_df = df[~df["Gene"].isin(genes_all_zero)].copy()
    after = len(filtered_df)
    rows_removed = before - after
    return filtered_df, genes_all_zero, rows_removed


def add_enrichment(
    agg_df: pd.DataFrame,
    gene_col: str = "Gene",
    value_col: str = "avg_nCPM",
    out_col: str = "Enrichment score",
    min_background: float = 1e-3,
    min_expression: float = 0.0,
    min_count: int = 2,
    pseudocount: Optional[float] = None,
    pseudocount_to_numerator: bool = False,
    clip_max: Optional[float] = None,
) -> pd.DataFrame:
    """Compute enrichment per row = current / mean(other cell types of the same gene)."""
    df = agg_df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    gene_sums = df.groupby(gene_col)[value_col].transform("sum")
    gene_counts = df.groupby(gene_col)[value_col].transform("count")

    denom_counts = gene_counts - 1
    avg_other = (gene_sums - df[value_col]) / denom_counts
    avg_other = avg_other.mask(denom_counts <= 0, np.nan)

    if pseudocount is not None:
        avg_other = avg_other + pseudocount

    numer = df[value_col]
    if pseudocount is not None and pseudocount_to_numerator:
        numer = numer + pseudocount

    denom = np.maximum(avg_other, min_background)
    numer = numer.where(numer >= min_expression, np.nan)

    df[out_col] = np.divide(
        numer, denom,
        out=np.full(df.shape[0], np.nan, dtype=float),
        where=(denom > 0)
    )
    df.loc[avg_other.isna(), out_col] = np.nan

    df[out_col] = df[out_col].where(gene_counts >= min_count, np.nan)

    if clip_max is not None:
        df[out_col] = df[out_col].clip(upper=clip_max)

    # Safe log2: add epsilon + mask zeros to NaN
    EPS = 1e-12
    log2_vals = np.log2(df[out_col].clip(lower=EPS))
    log2_vals = pd.Series(log2_vals, index=df.index).mask(df[out_col] <= 0)
    df["log2_enrichment"] = log2_vals
    return df


def batch_normalize_if_needed(df: pd.DataFrame, value_col: str, batch_col: Optional[str], batch_normalize: str = "none") -> pd.DataFrame:
    """Median-scale the value column per batch to align batch medians with the global median."""
    if not batch_col or batch_normalize == "none" or batch_col not in df.columns:
        return df

    out = df.copy()
    # Global median across all rows
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")
    global_median = out[value_col].median()
    if pd.isna(global_median) or global_median == 0:
        return out

    # Per-batch medians
    batch_medians = (
        out.groupby(batch_col, dropna=False)[value_col]
           .median()
           .rename("_batch_median")
    )

    # Map per-row scale = global_median / batch_median
    out = out.merge(batch_medians, on=batch_col, how="left")
    scale = np.where((out["_batch_median"].notna()) & (out["_batch_median"] != 0),
                     global_median / out["_batch_median"], 1.0)
    out[value_col] = out[value_col] * scale
    out.drop(columns=["_batch_median"], inplace=True)
    return out
def aggregate_with_celltype(
    df: pd.DataFrame,
    gene_col: str,
    gene_name_col: str,
    cell_type_col: str,
    value_col: str,
    cluster_col: Optional[str],
    weighted: bool,
    weight_col: Optional[str] = "Read count",
    cluster_aggregate: str = "mean",  # used when weighted is off: mean|median
) -> pd.DataFrame:
    """
    Aggregate to one row per (Gene × Gene name × Cell type) across clusters.
    If weighted=True and weight_col exists: weighted mean Σ(nCPM*w)/Σ(w).
    Else: unweighted mean or median across clusters.
    """
    df = df.copy()
    group_cols = [gene_col, gene_name_col, cell_type_col]

    if weighted and weight_col and (weight_col in df.columns):
        w = pd.to_numeric(df[weight_col], errors="coerce").fillna(1.0)
        vals = pd.to_numeric(df[value_col], errors="coerce")
        df["__prod__"] = vals * w
        agg = (
            df.groupby(group_cols, as_index=False)
              .agg(
                  avg_nCPM=("__prod__", "sum"),
                  weight_sum=(weight_col, "sum"),
                  clusters_used=(cluster_col, "nunique") if (cluster_col and cluster_col in df.columns) else (value_col, "count")
              )
        )
        agg["avg_nCPM"] = agg["avg_nCPM"] / agg["weight_sum"].replace(0, np.nan)
        df.drop(columns=["__prod__"], inplace=True)
    else:
        # Unweighted aggregation across clusters
        func = "mean" if cluster_aggregate == "mean" else "median"
        agg = (
            df.groupby(group_cols, as_index=False)
              .agg(
                  avg_nCPM=(value_col, func),
                  clusters_used=(cluster_col, "nunique") if (cluster_col and cluster_col in df.columns) else (value_col, "count")
              )
        )
    return agg

# ---------- Yanai's τ (specificity) ----------

def gene_specificity_tau(sub: pd.DataFrame, expr_col: str = "avg_nCPM") -> float:
    """Compute Yanai's τ for one gene from its avg_nCPM across cell types (0..1)."""
    vals = pd.to_numeric(sub[expr_col], errors="coerce").fillna(0.0).to_numpy()
    if len(vals) == 0:
        return np.nan
    m = vals.max()
    if m <= 0:
        return 0.0
    y = vals / m  # normalize by max
    K = len(vals)
    tau = (np.sum(1.0 - y)) / (K - 1) if K > 1 else 1.0
    return float(tau)


def gene_specificity_tau_series(s: pd.Series) -> float:
    """Yanai's τ from a Series of avg_nCPM values across cell types (0..1)."""
    vals = pd.to_numeric(s, errors="coerce").fillna(0.0).to_numpy()
    if len(vals) == 0:
        return np.nan
    m = vals.max()
    if m <= 0:
        return 0.0
    y = vals / m
    K = len(vals)
    tau = (np.sum(1.0 - y)) / (K - 1) if K > 1 else 1.0
    return float(tau)


def compute_tau(agg_df: pd.DataFrame, gene_col: str = "Gene", expr_col: str = "avg_nCPM") -> pd.DataFrame:
    """Compute τ per gene and merge as a report column specificity_tau."""
    df = agg_df.copy()
    tau_series = df.groupby(gene_col, group_keys=False)[expr_col].apply(gene_specificity_tau_series)
    tau_df = tau_series.rename("specificity_tau").reset_index()
    df = df.merge(tau_df, on=gene_col, how="left")
    return df


def apply_tau_filter(
    agg_df: pd.DataFrame,
    gene_col: str = "Gene",
    min_specificity: Optional[float] = None,
) -> Tuple[pd.DataFrame, int, int]:
    """
    Drop genes with τ < min_specificity. Assumes specificity_tau is present.
    Returns (filtered_df, n_total_genes, n_dropped_genes).
    """
    if min_specificity is None:
        total = agg_df[gene_col].nunique()
        return agg_df, total, 0
    df = agg_df.copy()
    tau_per_gene = df[[gene_col, "specificity_tau"]].drop_duplicates(subset=[gene_col])
    n_total = tau_per_gene[gene_col].nunique()
    keep_genes = tau_per_gene.loc[tau_per_gene["specificity_tau"] >= min_specificity, gene_col]
    filtered_df = df[df[gene_col].isin(keep_genes)].copy()
    n_dropped = n_total - keep_genes.nunique()
    return filtered_df, n_total, n_dropped

# ---------- Top N helpers ----------

def top_n_overall(df: pd.DataFrame, sort_by: str, n: int) -> pd.DataFrame:
    return df.sort_values(by=sort_by, ascending=False).head(n).copy()


def top_n_per_cell_type(df: pd.DataFrame, cell_type_col: str, sort_by: str, n: int) -> pd.DataFrame:
    out_frames = []
    for ct, sub in df.groupby(cell_type_col):
        sub_sorted = sub.sort_values(by=sort_by, ascending=False).head(n).copy()
        sub_sorted["rank_in_cell_type"] = range(1, len(sub_sorted) + 1)
        out_frames.append(sub_sorted)
    if out_frames:
        return pd.concat(out_frames, axis=0, ignore_index=True)
    else:
        return pd.DataFrame(columns=list(df.columns) + ["rank_in_cell_type"])

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description="Cell-type–aware enrichment with optional weighting and τ report/filter/penalize — v1.4."
    )

    # I/O
    parser.add_argument("--input-file", type=str, default="rna_single_cell_cluster.tsv",
                        help="Path to input TSV containing single-cell cluster data.")
    parser.add_argument("--output-file", type=str, default="celltype_enrichment.tsv",
                        help="Path for the full enrichment output TSV.")
    parser.add_argument("--top-n", type=int, default=100,
                        help="Top N rows to save overall (default: 100).")
    parser.add_argument("--per-cell-type-top-n", type=int, default=20,
                        help="Top N per cell type to export (0 disables).")

    # Column names (defaults match your example)
    parser.add_argument("--gene-col", type=str, default="Gene")
    parser.add_argument("--gene-name-col", type=str, default="Gene name")
    parser.add_argument("--cell-type-col", type=str, default="Cell type")
    parser.add_argument("--cluster-col", type=str, default="Cluster")
    parser.add_argument("--batch-col", type=str, default="Cell type",
                        help="Batch column used for median scaling (default: Cell type).")
    parser.add_argument("--value-col", type=str, default="nCPM")
    parser.add_argument("--weight-col", type=str, default="Read count",
                        help="Weight column used for weighted aggregation (default: Read count).")

    # Weighting options
    parser.add_argument("--weighted", type=str, choices=["on", "off"], default="on",
                        help="Use weighted aggregation across clusters (on/off).")
    parser.add_argument("--cluster-aggregate", type=str, choices=["mean", "median"], default="mean",
                        help="When --weighted off, choose mean or median across clusters.")

    # Filters and safeguards
    parser.add_argument("--min-clusters", type=int, default=None,
                        help="Keep only rows where clusters_used >= this integer.")
    parser.add_argument("--treat-nan-as-zero", action="store_true",
                        help="Treat NaN as zero when deciding 'no expression' genes to drop.")
    parser.add_argument("--min-expr-threshold", type=float, default=0.0,
                        help="Filter rows with avg_nCPM < threshold before enrichment (default: 0).")

    # Enrichment parameters
    parser.add_argument("--min-background", type=float, default=1e-3,
                        help="Minimum denominator for enrichment (default: 1e-3).")
    parser.add_argument("--min-expression", type=float, default=0.0,
                        help="Minimum numerator expression to compute enrichment (default: 0).")
    parser.add_argument("--min-count", type=int, default=2,
                        help="Require >= this many entries per gene for enrichment.")
    parser.add_argument("--pseudocount", type=float, default=None,
                        help="Optional pseudocount added to denominator; set a small value like 0.01.")
    parser.add_argument("--pseudocount-to-numerator", action="store_true",
                        help="Also add pseudocount to numerator.")
    parser.add_argument("--clip-max", type=float, default=None,
                        help="Optional cap on enrichment score (e.g., 100).")
    parser.add_argument("--sort-by", type=str, choices=[
        "Enrichment score",
        "log2_enrichment",
        "Enrichment score (tau penalized)",
        "log2_enrichment_penalized"
    ], default="log2_enrichment",
        help="Column used for sorting outputs.")

    # Batch normalization
    parser.add_argument("--batch-normalize", type=str, choices=["none", "median_scale"],
                        default="median_scale", help="Per-batch normalization method.")

    # τ report/filter/penalize
    parser.add_argument("--specificity-mode", type=str, choices=["off", "filter", "penalize"], default="off",
                        help="How to use Yanai's τ: off, filter (drop genes below threshold), or penalize enrichment.")
    parser.add_argument("--min-specificity", type=float, default=None,
                        help="Threshold for τ (0..1). Used by --specificity-mode filter|penalize. Example: 0.8.")

    args = parser.parse_args()

    # -------- Load input --------
    print("\033[33mLoading input...\033[0m")
    df = pd.read_csv(args.input_file, sep="\t")

    # Sanity check required columns
    required = [args.gene_col, args.gene_name_col, args.cell_type_col, args.value_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Input is missing required columns: {missing}")

    # -------- Coerce numeric on value_col and report non-numeric --------
    df[args.value_col], non_numeric = coerce_numeric(df[args.value_col])
    if non_numeric:
        print(f"\033[31mNon-numeric values in {args.value_col}: {non_numeric}\033[0m")
    else:
        print("\033[32mGOOD DATA nCPM. NO Non numeric\033[0m")

    # -------- Optional batch normalization (per cell type by default) --------
    if args.batch_col and args.batch_col in df.columns:
        print(f"\033[33mApplying batch normalization per '{args.batch_col}'...\033[0m")
        df = batch_normalize_if_needed(
            df, value_col=args.value_col,
            batch_col=args.batch_col,
            batch_normalize=args.batch_normalize
        )
    else:
        if args.batch_col:
            print(f"\033[33mBatch column '{args.batch_col}' not found; skipping normalization.\033[0m")

    # -------- Aggregation across clusters within each (Gene × Cell type) --------
    print("\033[33mAggregating across clusters within each (Gene × Cell type)...\033[0m")
    agg_df = aggregate_with_celltype(
        df,
        gene_col=args.gene_col,
        gene_name_col=args.gene_name_col,
        cell_type_col=args.cell_type_col,
        value_col=args.value_col,
        cluster_col=args.cluster_col if args.cluster_col in df.columns else None,
        weighted=(args.weighted == "on"),
        weight_col=args.weight_col if args.weight_col in df.columns else None,
        cluster_aggregate=args.cluster_aggregate
    )

    # -------- Optional filter by min clusters --------
    if args.min_clusters is not None:
        before_n = len(agg_df)
        agg_df = agg_df[agg_df["clusters_used"] >= args.min_clusters].copy()
        after_n = len(agg_df)
        print(
            f"\033[33mFiltered by clusters_used >= {args.min_clusters}. "
            f"Dropped {before_n - after_n} row(s); {after_n} row(s) remain.\033[0m"
        )

    # -------- Pre-enrichment expression threshold --------
    if args.min_expr_threshold > 0:
        before_n = len(agg_df)
        agg_df = agg_df[agg_df["avg_nCPM"] >= args.min_expr_threshold].copy()
        after_n = len(agg_df)
        print(
            f"\033[33mFiltered rows with avg_nCPM < {args.min_expr_threshold}. "
            f"Dropped {before_n - after_n} row(s); {after_n} row(s) remain.\033[0m"
        )

    # -------- Drop genes with no expression across all cell types --------
    print("\033[33mDropping genes with no expression across all cell types...\033[0m")
    agg_df, dropped_genes, rows_removed = drop_genes_with_no_expression(
        agg_df, expr_col="avg_nCPM", treat_nan_as_zero=args.treat_nan_as_zero
    )
    print(
        f"\033[31mDropped {len(dropped_genes)} gene(s), removing {rows_removed} row(s).\033[0m"
    )

    # -------- Compute τ and merge as a report column --------
    print("\033[33mComputing Yanai's τ per gene...\033[0m")
    agg_df = compute_tau(agg_df, gene_col=args.gene_col, expr_col="avg_nCPM")

    # -------- Specificity mode: filter or penalize --------
    if args.specificity_mode == "filter" and args.min_specificity is not None:
        print(f"\033[33mApplying τ filtering (threshold={args.min_specificity})...\033[0m")
        agg_df, n_total, n_dropped = apply_tau_filter(
            agg_df, gene_col=args.gene_col, min_specificity=args.min_specificity
        )
        print(
            f"\033[33mτ filter: {n_dropped} gene(s) dropped out of {n_total}. "
            f"{agg_df[args.gene_col].nunique()} gene(s) remain.\033[0m"
        )
    elif args.specificity_mode == "filter" and args.min_specificity is None:
        print("\033[33mWARNING: --specificity-mode filter set but --min-specificity not provided; skipping filter.\033[0m")

    # -------- Unique cell types list --------
    n_cell_types = agg_df[args.cell_type_col].dropna().nunique()
    unique_cell_types = sorted(agg_df[args.cell_type_col].dropna().unique().tolist())
    pd.Series(unique_cell_types, name=args.cell_type_col).to_csv("unique_cell_types.tsv", sep="\t", index=False)
    print(f"\033[32mNumber of unique cell types: {n_cell_types}\033[0m")

    # -------- Count unique genes --------
    n_unique_genes = agg_df[args.gene_col].dropna().nunique()
    print(f"\033[32mNumber of unique genes remaining: {n_unique_genes}\033[0m")

    # -------- Enrichment --------
    print("\033[33mCalculating Enrichment Scores...\033[0m")
    agg_df = add_enrichment(
        agg_df=agg_df,
        gene_col=args.gene_col,
        value_col="avg_nCPM",
        out_col="Enrichment score",
        min_background=args.min_background,
        min_expression=args.min_expression,
        min_count=args.min_count,
        pseudocount=args.pseudocount,
        pseudocount_to_numerator=args.pseudocount_to_numerator,
        clip_max=args.clip_max
    )
    print("\033[33mDone calculating.\033[0m")

    # -------- Penalize enrichment by τ (optional) --------
    if args.specificity_mode == "penalize":
        print("\033[33mApplying τ penalization...\033[0m")
        penalty = agg_df["specificity_tau"].clip(lower=0.0, upper=1.0)
        agg_df["Enrichment score (tau penalized)"] = agg_df["Enrichment score"] * penalty
        EPS = 1e-12
        log2p = np.log2(agg_df["Enrichment score (tau penalized)"].clip(lower=EPS))
        log2p = pd.Series(log2p, index=agg_df.index).mask(agg_df["Enrichment score (tau penalized)"] <= 0)
        agg_df["log2_enrichment_penalized"] = log2p

    # -------- Single cell-type genes --------
    gene_celltype_counts = agg_df.groupby(args.gene_col)[args.cell_type_col].transform("nunique")
    agg_df["single_cell_type_gene"] = (gene_celltype_counts == 1)
    min_ct_per_gene = int(gene_celltype_counts.min()) if len(gene_celltype_counts) else 0
    print(
        f"\033[33mMinimum number of cell types per gene (across genes): {min_ct_per_gene}\033[0m"
    )

    single_cell_rows = agg_df[agg_df["single_cell_type_gene"].fillna(False)].copy()
    if not single_cell_rows.empty:
        n_genes_single = single_cell_rows[args.gene_col].nunique()
        print(f"\033[32mFound {n_genes_single} gene(s) expressed in exactly one cell type- according to raw data.\033[0m")
        single_cell_rows.to_csv("single_cell_type_gene_rows.tsv", sep="\t", index=False)
    else:
        print("\033[33mNo genes were found that are only expressed in one cell type.\033[0m")

    # -------- Sort and save full --------
    sort_col = args.sort_by
    agg_df_sorted = agg_df.sort_values(by=sort_col, ascending=False)
    agg_df_sorted.to_csv(args.output_file, sep="\t", index=False)
    print(f"\033[33mSaved full table: {args.output_file}\033[0m")

    # -------- Top-N overall --------
    print("\033[33mSaving top-N overall...\033[0m")
    top_overall = top_n_overall(agg_df_sorted, sort_by=sort_col, n=args.top_n)
    suffix = 'log2' if 'log2' in sort_col else ('penalized' if 'penalized' in sort_col else 'enrichment')
    top_overall_file = f"top{args.top_n}_{suffix}.tsv"
    top_overall.to_csv(top_overall_file, sep="\t", index=False)
    print(f"\033[33mSaved: {top_overall_file}\033[0m")

    # -------- Top-N per cell type (skip if N=0) --------
    if args.per_cell_type_top_n and args.per_cell_type_top_n > 0:
        print("\033[33mSaving top-N per cell type...\033[0m")
        top_pct = top_n_per_cell_type(agg_df, cell_type_col=args.cell_type_col, sort_by=sort_col, n=args.per_cell_type_top_n)
        top_pct_suffix = 'log2' if 'log2' in sort_col else ('penalized' if 'penalized' in sort_col else 'enrichment')
        top_pct_file = f"top_per_cell_type_{args.per_cell_type_top_n}_{top_pct_suffix}.tsv"
        top_pct.to_csv(top_pct_file, sep="\t", index=False)
        print(f"\033[33mSaved: {top_pct_file}\033[0m")
    else:
        print("\033[33mPer-cell-type top-N export disabled (N=0).\033[0m")


if __name__ == "__main__":
    main()
```
## 3-III Celltype enrichment runner
I ran it with the following runner

run_celltype_enrichment_v1_4.sh
```
#!/usr/bin/env bash
# run_celltype_enrichment_with_options.sh (v1.4 with τ report & penalization)
# Usage:
#   ./run_celltype_enrichment_with_options.sh [overrides...]
# Examples:
#   ./run_celltype_enrichment_with_options.sh --specificity-mode filter --min-specificity 0.8
#   ./run_celltype_enrichment_with_options.sh --specificity-mode penalize --min-specificity 0.8 --sort-by "log2_enrichment_penalized"

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

args=(
  # Path to your input TSV containing: Gene, Gene name, Cell type, Cluster, Read count, nCPM
  --input-file integrated_filtered.tsv
  # Path for the full enrichment output table (TSV)
  --output-file celltype_enrichment.tsv
  # Number of top rows to export overall (after sorting)
  --top-n 100
  # Number of top rows to export per cell type (set 0 to disable)
  --per-cell-type-top-n 20

  # Column name for gene ID (e.g., ENSG...)
  --gene-col "Gene"
  # Column name for gene symbol/display name
  --gene-name-col "Gene name"
  # Column name for cell type labels
  --cell-type-col "Cell type"
  # Column name for cluster IDs (replicates within a cell type)
  --cluster-col "Cluster"
  # Batch column used for median scaling (default: per cell type)
  --batch-col "Cell type"
  # Column name for expression values to aggregate
  --value-col "nCPM"
  # Column name for weights used in weighted aggregation
  --weight-col "Cell count"

  # Toggle weighted aggregation across clusters: on=weighted mean, off=unweighted
  --weighted on
  # If weighted is off, choose how to aggregate clusters: mean or median
  #--cluster-aggregate mean

  # Minimum number of clusters required to keep a (Gene × Cell type) row
  --min-clusters 2
  # Drop rows with avg_nCPM below this threshold before enrichment. This helps you to deal with weird values generated by 0 expression 
  --min-expr-threshold 0.00
  # Treat NaN as zero when deciding to drop genes with no expression (comment to disable)
  # --treat-nan-as-zero

  # Floor for the denominator in enrichment to avoid tiny values
  --min-background 1e-3
  # Minimum numerator expression required to compute enrichment
  --min-expression 0.001
  # Minimum number of cell-type entries per gene to compute enrichment
  --min-count 2
  # Add a small constant to stabilize denominators (uncomment to enable)
  --pseudocount 0.001
  # Also add pseudocount to numerator (use with --pseudocount)
  # --pseudocount-to-numerator
  # Cap extremely large enrichment ratios (uncomment to enable)
  # --clip-max 100

  # Sort outputs by: raw/log2 or penalized variants
  --sort-by "log2_enrichment_penalized"

  # Apply median scaling per batch (here: per cell type)
  --batch-normalize median_scale

  # Yanai's τ usage:
  #   off      → only report `specificity_tau` (no filtering/penalization)
  #   filter   → drop genes whose τ < threshold (strict specificity)
  #   penalize → multiply enrichment by τ (keeps genes but down-ranks broad ones)
  --specificity-mode penalize
  # τ threshold (0..1) used by filter/penalize modes; e.g., 0.8 for high specificity
  --min-specificity 0.8
)

# Allow overrides on the CLI (last wins)
args+=("$@")

# Call the Python script
python3 "${script_dir}/celltype_enrichment_v1_4.py" "${args[@]}"

```
## 3-IV Run command
Then I ran it like
```
./run_celltype_enrichment_v1_4.sh --input-file combined_expression_data_filtered.tsv --output-file enrichment_values_for_filtered_celltypes.tsv --min-clusters 1 --min-count 50 --specificity-mode penalize --min-specificity 1
```
# 4) Visualizing Enrichment Distributions

## 4-I Introduction

This section provides a simple CLI tool to visualize the distribution of values in any numeric column (e.g., log2_enrichment_penalized) from a TSV/CSV file. It adds:

A horizontal line at 0
Grey shading for the depletion (negative) zone
Clear labels for Enrichment and Depletion
The ability to highlight near-zero rows with red vertical markers
A function to report the % of positive values
Optional export of the near-zero rows to a CSV

Why this matters

log₂ > 0 → Enriched (observed higher than expected)
log₂ < 0 → Depleted (observed lower than expected)
log₂ = 0 → Neutral (observed ≈ expected; enrichment score = 1)

This script makes it easy to see that boundary and inspect data near the neutral point.

## 4-II Script

plot_distribution.py
```py

#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_positive_percentage(series):
    total = len(series)
    positive_count = (series > 0).sum()
    return (positive_count / total) * 100

def get_rows_near_zero(df, column, n):
    # Drop NA values and sort by absolute distance from zero
    df_clean = df.loc[df[column].notna()].copy()
    df_clean["abs_val"] = df_clean[column].abs()
    return df_clean.sort_values("abs_val").drop(columns="abs_val").head(n)

def main():
    parser = argparse.ArgumentParser(description="Visualize distribution and analyze enrichment values.")
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV/TSV file")
    parser.add_argument("-c", "--column", required=True, help="Column name to plot distribution")
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter (default: tab)")
    parser.add_argument("-o", "--output", default=None, help="Output plot filename")
    parser.add_argument("--near-zero", type=int, default=0, help="Show N rows closest to zero")
    parser.add_argument("--export-near-zero", default=None, help="Export near-zero rows to CSV file")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input, delimiter=args.delimiter)

    if args.column not in df.columns:
        print(f"Error: Column '{args.column}' not found. Available columns: {list(df.columns)}")
        return

    # Positive percentage
    pct_positive = calculate_positive_percentage(df[args.column])
    print(f"Percentage of rows with positive {args.column}: {pct_positive:.2f}%")

    # Plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 7))
    ax = sns.histplot(df[args.column], kde=True, color="blue", bins=30)

    # Vertical line at 0
    plt.axvline(0, color="black", linestyle="--", linewidth=1)

    # Get axis limits after plotting
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    # Shade depletion zone
    plt.axvspan(xmin, 0, color="grey", alpha=0.2)

    # Labels for zones
    plt.text(xmin + (xmax - xmin) * 0.05, ymax * 0.92, "Depletion", color="grey", fontsize=12)
    plt.text(xmin + (xmax - xmin) * 0.70, ymax * 0.92, "Enrichment", color="blue", fontsize=12)

    # Highlight near-zero rows if requested
    if args.near_zero > 0:
        near_zero_rows = get_rows_near_zero(df, args.column, args.near_zero)
        for val in near_zero_rows[args.column]:
            plt.axvline(val, color="red", linestyle=":", linewidth=1)
        print(f"\nTop {args.near_zero} rows closest to zero in {args.column}:")
        print(near_zero_rows)

        # Export if requested
        if args.export_near_zero:
            near_zero_rows.to_csv(args.export_near_zero, index=False)
            print(f"Near-zero rows exported to {args.export_near_zero}")

    # Title and labels
    plt.title(f"Distribution of {args.column}")
    plt.xlabel(args.column)
    plt.ylabel("Frequency")

    # Save plot
    output_file = args.output if args.output else f"{args.column}_distribution.png"
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved as {output_file}")

if __name__ == "__main__":
    main()
```
## 4-III CLI Help
```txt

usage: plot_distribution.py [-h] -i INPUT -c COLUMN [-d DELIMITER] [-o OUTPUT]
                            [--near-zero NEAR_ZERO] [--export-near-zero EXPORT_NEAR_ZERO]

Visualize distribution and analyze enrichment values.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input CSV/TSV file
  -c COLUMN, --column COLUMN
                        Column name to plot distribution
  -d DELIMITER, --delimiter DELIMITER
                        Delimiter (default: tab)
  -o OUTPUT, --output OUTPUT
                        Output plot filename
  --near-zero NEAR_ZERO
                        Show N rows closest to zero
  --export-near-zero EXPORT_NEAR_ZERO
                        Export near-zero rows to CSV file
```
## 4-IV Run

#* I ran it like following
```bash
python plot_distribution.py -i enrichment_values_for_filtered_celltypes.tsv -c log2_enrichment_penalized --near-zero 10 -o log2_enrichment_penalized_distribution.png
```
This code plots the distribution, Give you the top percentage of rows with positive values and shows you the rows around that value
```bash
Percentage of rows with positive log2_enrichment_penalized: 18.17%
```
With this data, in the ranking step below (section 5), I am 

1. dropping all the negative values and 0s as they represent cell type*gene combinations with equal or less than expression to the background
2. Checking the presense of expression of genes in cell types in top x% of the enrichment values FROM THE REMAINING values so I only get gene expression which is significantly higher than background 

log2_enrichment_penalized_distribution.png

<img width="1000" height="700" alt="log2_enrichment_penalized_distribution" src="https://github.com/user-attachments/assets/ad770994-6ef6-4b93-a0b4-f0fe4915e96b" />


# 5) Add celltype info layers

## 5-I Introduction

My enrichment value table only has cell type information. But it can be useful to have information like cell type group and cell type class. I am adding those layers to the dataset here

HPA has these information on the table here:
```url
https://www.proteinatlas.org/download/tsv/rna_single_cell_type_cell_types.tsv.zip
```
I can use the same merging script from section 1

## 5-II Run command

I Run it like following

```bash
python ../1.Combining_datasets/1.Raw_data/merge_tsv_by_keys.py \
  --left rna_single_cell_type_cell_types_with_added_groups.txt \
  --right rna_single_cell_type_cell_types.tsv \
  --left-keys "Cell type" \
  --right-keys "Cell type" \
  --right-cols "Cell type group,Cell type class" \
  --out  enrich_values_with_cell_class_data.tsv
```
## 5-III Fixing celltypes with missing groups

As I am changing filtering values as needed, sometimes 'rna_single_cell_type_cell_types.tsv' does not have all the 'cell type groups' and 'cell type classes' for all the cell types in enrichment table. Therefore I have to look for rows with missing values for 'Cell type group' and ' Cell type class' with following command
```bash
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Cell type group") col=i; next} $col=="" {print}' enrich_values_with_cell_class_data.tsv | cut -f 3 | sort -u
```
Then I added information of missing groups and classes to a new file named 'rna_single_cell_type_cell_types_with_added_groups.txt' and re mapped cell type groups and classses to rna_single_cell_type_cell_types_with_added_groups.txt to get a complete set.

Following are the rows I added
```tsv
cycling epithelial cells	proliferating cells	stem and proliferating cells
cycling immune cells	proliferating cells	blood and immune cells
cycling salivary epithelial cells	proliferating cells	specialized epithelial cells
fallopian basal stem cells	stem cells	stem and proliferating cells
fetal immune cells	other immune cells	blood and immune cells
glandular and luminal cells	glandular and luminal cells	glandular epithelial cells
intestinal stem cells	stem cells	stem and proliferating cells
mesangial cells	mural cells	endothelial and mural cells
multipotent progenitors	stem cells	stem and proliferating cells
mural & epithelial cells	mural cells	endothelial and mural cells
mural cells	mural cells	endothelial and mural cells
stromal cells	stromal cells	mesenchymal cells
```

Then ran it again with the updtaed file


```
python ../1.Combining_datasets/1.Raw_data/merge_tsv_by_keys.py \
  --left enrichment_values_for_filtered_celltypes.tsv \
  --right rna_single_cell_type_cell_types_with_added_groups.txt \
  --left-keys "Cell type" \
  --right-keys "Cell type" \
  --right-cols "Cell type group,Cell type class" \
  --out  enrich_values_with_cell_class_data.tsv
```

# 6) Rank genes on cell specific expresion

## 6-I Introduction
rank_genes.py ranks genes by enrichment and summarizes how many groups (e.g., Cell type, Cell type group, Cell type class) each gene appears in within a global top subset.
You have full control over:

Grouping column (via --presence-col) used for ranking and “rank within group”.
Uniqueness column (via --unique-col) used to count distinct values per gene.
Selection & sorting metrics (--top-col, --sorting-col) and top‑% subset size.
Output column names, filtering, and NA/negative handling.

This lets you switch between cell‑type specific, group‑level, or class‑level views without changing the code—just alter the CLI flags.

## 6-II Rank genes and estimate celltype counts script

rank_genes.py
```py

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


```
## 6-III CLI help
```txt


usage: rank_genes.py [-h] --input INPUT --output OUTPUT
                     [--top-percent TOP_PERCENT] [--min-top-rows MIN_TOP_ROWS]
                     [--top-col TOP_COL] [--sorting-col SORTING_COL]
                     [--gene-col GENE_COL] [--celltype-col CELLTYPE_COL]
                     [--presence-col PRESENCE_COL] [--unique]
                     [--unique-col UNIQUE_COL] [--drop-na] [--drop-negatives]
                     [--output-count-col OUTPUT_COUNT_COL]
                     [--output-list-col OUTPUT_LIST_COL]
                     [--output-overall-rank-col OUTPUT_OVERALL_RANK_COL]
                     [--output-rank-within-col OUTPUT_RANK_WITHIN_COL]
                     [--include-cols INCLUDE_COLS [INCLUDE_COLS ...]] [--verbose]

Rank genes globally by enrichment and count groups.

required arguments:
  --input INPUT
        Input TSV file (tab-delimited) containing columns like:
        Gene, Gene name, Cell type, log2_enrichment_penalized, etc.

  --output OUTPUT
        Output TSV file to write ranked results.

top subset selection:
  --top-percent TOP_PERCENT
        Percentage (0–100) used to form the global top subset (default: 25.0).

  --min-top-rows MIN_TOP_ROWS
        Minimum number of rows to retain in the top subset (default: 1).

  --top-col TOP_COL
        Column used to determine eligibility for the top subset
        (numeric; default: log2_enrichment_penalized).

  --sorting-col SORTING_COL
        Column used to sort descending before building the subset and
        as the secondary key in the final ranking (numeric; default: log2_enrichment_penalized).

entity & grouping:
  --gene-col GENE_COL
        Gene identifier column (default: Gene).

  --celltype-col CELLTYPE_COL
        Deprecated fallback for presence; used only if --presence-col is not set (default: "Cell type").

  --presence-col PRESENCE_COL
        Column defining the grouping used for "rank within group" and selection of rows
        (e.g., "Cell type", "Cell type group", "Cell type class"). If omitted, falls back to --celltype-col.

uniqueness controls:
  --unique
        Enable deduplication when counting groups per gene, using --unique-col.

  --unique-col UNIQUE_COL
        Column used to count unique values per gene (e.g., "Cell type group").
        If omitted, defaults to the value of --presence-col.

filters:
  --drop-na
        Drop rows where --top-col or --sorting-col are NaN after numeric coercion.

  --drop-negatives
        Drop rows where --top-col or --sorting-col are negative.

output controls:
  --output-count-col OUTPUT_COUNT_COL
        Name of the per-gene group count column (default: group_count).

  --output-list-col OUTPUT_LIST_COL
        Name of the per-gene group list column (default: group_list).

  --output-overall-rank-col OUTPUT_OVERALL_RANK_COL
        Name of the overall rank column (default: overall_rank).

  --output-rank-within-col OUTPUT_RANK_WITHIN_COL
        Name of the rank-within-group column (default: rank_within_group).

  --include-cols INCLUDE_COLS [INCLUDE_COLS ...]
        Extra columns to keep in the output (default includes:
        "Gene", "Gene name", "Cell type", "avg_nCPM", "specificity_tau",
        "Enrichment score (tau penalized)", "log2_enrichment_penalized").

misc:
  --verbose
        Print summary information to stdout.

  -h, --help
        Show this help message and exit.

```
## 6-IV Run command
#*I ran it like following*

#*I Ranked these first by 'Cell type'*

#*And only checked the prescense in top 95% of the enrichment values as lowest 5% could be close to background noise*

```
python rank_genes.py \
  --input enrich_values_with_cell_class_data.tsv --drop-na --drop-negatives\
  --output ranked_genes_unique_celltypes.tsv \
  --presence-col "Cell type" \
  --unique-col "Cell type" \
  --unique \
  --top-percent 95 \
  --top-col log2_enrichment_penalized \
  --sorting-col log2_enrichment_penalized \
  --gene-col Gene \
  --output-count-col top_percent_Cell_type_count \
  --output-list-col top_percent_Cell_types \
  --output-rank-within-col rank_within_Cell_type \
  --output-overall-rank-col overall_rank_by_Cell_type \
  --verbose

```

#*Then I rankd them again by the PRESCENCE IN 'Cell type group'*

This does no change the previous ranking. This just adds another column to it with ranks based on cell group instead cell type

STILL I DID NOT CALCULATE ENRICHMENT VALUES SEPERATELY FOR 'Cell type group'. ENRICHMENT IS STILL BASED ON ALL THE 'Cell Type' values

AND I AM SETTING --top-percent TO 100% HERE AS IT WAS ALREADY FILTERED IN THE INPUT FILE USED HERE. Same with --drop-na and negatives etc

```
python rank_genes.py \
  --input ranked_genes_unique_celltypes.tsv \
  --output ranked_genes_unique_celltypes_and_groups.tsv \
  --presence-col "Cell type group" \
  --unique-col "Cell type group" \
  --unique \
  --top-percent 100 \
  --top-col log2_enrichment_penalized \
  --sorting-col log2_enrichment_penalized \
  --gene-col Gene \
  --output-count-col top_percent_Cell_type_group_count \
  --output-list-col top_percent_Cell_type_groups \
  --output-rank-within-col rank_within_Cell_type_group \
  --output-overall-rank-col overall_rank_by_Cell_type_group \
  --verbose
```
#*Then I rankd them again by the PRESCENCE IN 'Cell type class'*
```
python rank_genes.py \
  --input ranked_genes_unique_celltypes_and_groups.tsv \
  --output ranked_genes_unique_celltypes_and_groups_and_classes.tsv \
  --presence-col "Cell type class" \
  --unique-col "Cell type class" \
  --unique \
  --top-percent 100 \
  --top-col log2_enrichment_penalized \
  --sorting-col log2_enrichment_penalized \
  --gene-col Gene \
  --output-count-col top_percent_Cell_type_class_count \
  --output-list-col top_percent_Cell_type_classes \
  --output-rank-within-col rank_within_Cell_type_class \
  --output-overall-rank-col overall_rank_by_Cell_type_class \
  --verbose
```
# 7) Cluster data integration

## 7-I Introduction

#* First I added cluster categories (to add as a filter later on section 8) with a minimum limit with this script. You can select a custom number of clusters to filter by. Here I use 

1 or below and

2 or above as limits

## 7_II Cluster categories script

categorize_column.py
```py

import pandas as pd
import sys

def main():
    if len(sys.argv) != 6:
        print("Usage: python categorize_column.py <input_file> <output_file> <selected_column> <selected_number> <result_column>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    selected_column = sys.argv[3]
    selected_number = float(sys.argv[4])
    result_column = sys.argv[5]

    # Load input file (auto-detect delimiter)
    try:
        df = pd.read_csv(input_file, sep=None, engine="python")
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    if selected_column not in df.columns:
        print(f"Error: column '{selected_column}' not found.")
        sys.exit(1)

    # Create categorized column
    df[result_column] = df[selected_column].apply(
        lambda x: "1 or below" if float(x) < selected_number else "2 or higher"
    )

    # Save output *as TSV*
    try:
        df.to_csv(output_file, sep="\t", index=False)
        print(f"TSV output written to {output_file}")
    except Exception as e:
        print(f"Error writing output: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```
## 7-III CLI help
```
python categorize_column.py <input_file> <output_file> <selected_column> <selected_number> <result_column>
```
## 7-IV Run command

I ran it like following
```bash
python categorize_column.py ranked_genes_unique_celltypes_and_groups_and_classes.tsv ranked_genes_with_cluster_categories.tsv c
lusters_used 2 cluster_limit
```

# 8) Protein type data integration

## 8-I Introduction

I found data on protein types from HPA:

"The human secretome comprises all proteins that are potentially secreted from the cell. These proteins were identified based on UniProt annotations of subcellular locations, along with computational predictions of signal peptides and transmembrane regions. The signal peptide is found in most secreted proteins, but also appears in some classes of membrane proteins. Therefore, the presence of transmembrane regions can be used to distinguish membrane proteins from secreted proteins. A whole-proteome scan of all Ensembl transcripts was performed using majority decision methods for signal peptide prediction (MDSEC) and membrane region prediction (MDM). Proteins with a predicted signal peptide by MDSEC and no predicted transmembrane region by MDM were considered secreted. Using this combined approach, the analysis predicts that 2520 genes (12% of all human protein-coding genes) encode at least one predicted secreted transcript or have at least one isoform with the subcellular location "Secreted" in UniProt." 

HPA: https://www.proteinatlas.org/humanproteome/tissue/secretome#classification_of_the_secretome

And I downloaded the dataset from
```url
https://www.proteinatlas.org/search/sa_location%3ASecreted+-+unknown+location%2CSecreted+in+brain%2CSecreted+in+female+reproductive+system%2CSecreted+in+male+reproductive+system%2CSecreted+in+other+tissues%2CSecreted+to+blood%2CSecreted+to+digestive+system%2CSecreted+to+extracellular+matrix%2CIntracellular+and+membrane%2CImmunoglobulin+genes
```
It had protein type data in the column 'Protein type' as the first comma seperated value that starts with 'Predicted'. I extracted those preotei type data with other needed columns with the following script

## 8-II Extract essential data

### Script

extract_regex.py
```py

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
```

### CLI help

```txt
USAGE
-----
    python extract_regex.py [options]

OPTIONS
-------

Input / Output
--------------
--input <path>           (Required) Path to input TSV/CSV file.
--output <path>          (Required) Path to output file.
--file-delim <char>      File delimiter. Default: \t
--output-delim <char>    Output file delimiter. Default: \t
--encoding <encoding>    File encoding. Default: utf‑8.

Column Selection
----------------
--id-cols <col1 col2...> (Required) Columns to keep in output (e.g., Gene Ensembl).
--target-col <column>    (Required) Column containing the delimited values to parse.

Token Parsing
-------------
--value-sep <char>       Separator inside the target column. Default: ,
--pattern <regex>        (Required) Regex pattern to match tokens.
                         Examples:
                           "Predicted|Secreted"
                           "^(Predicted|Secreted)"
                           ".*(Predicted|Secreted).*"
--nth <int>              Select the n‑th matching token (1‑based). Default: 1.
--all-matches			 Gets all matches in addition to n-th
--ignore-case            Enable case‑insensitive regex matching.
--match-anywhere         Use re.search() instead of re.match().
                         Default behavior uses re.match() (anchored at start).

Performance
-----------
--chunksize <int>        Process file in chunks (e.g., 200000) for large datasets.

BEHAVIOR SUMMARY
----------------
• The target column is split into tokens using the delimiter (default comma).
• Tokens are trimmed.
• Regex is applied to each token:
    - Default: re.match() → matches only at the start of the token.
    - With --match-anywhere: re.search() → matches anywhere in token.
• The script extracts the n‑th match (1‑based index).
• Output contains:
    - All selected --id-cols
    - A new column named:
        <target-col>__regex_nth<n>

EXAMPLES
--------

1) First token starting with “Predicted” or “Secreted”
------------------------------------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out.tsv \
        --file-delim $'\t' \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern 'Predicted|Secreted' \
        --nth 1 \
        --ignore-case

2) Second matching token
------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out2.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --pattern 'Predicted|Secreted' \
        --nth 2 \
        --ignore-case

3) Match anywhere, not just start
---------------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out_anywhere.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern '(Predicted|Secreted)' \
        --match-anywhere \
        --ignore-case

4) Large files: chunked processing
----------------------------------
    python extract_regex.py \
        --input bigfile.tsv \
        --output out_big.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --pattern 'Predicted|Secreted' \
        --chunksize 250000

NOTES
-----
• Regex patterns must be quoted properly in shell environments.
• For TSV files, prefer --file-delim $'\t' to ensure correct parsing.
• Missing or empty target values produce empty output tokens.
```
### Run
```bash
python extract_regex.py \
  --input secretome_dataset.tsv \
  --output genes_with_protein_type.tsv \
  --file-delim $'\t' \
  --id-cols Gene Ensembl \
  --target-col "Protein class" \
  --value-sep "," \
  --pattern 'Predicted' \
  --ignore-case \
  --all-matches \
  --join-matches-sep ": "
```
## 8-III Merging protein type data with enrichment

I Started by renaming the columns of genes_with_protein_type.tsv in a meaningful way to make it easier to map

First I 

### made a mapping tsv 

file for column names with new column names

new_header.txt
```txt
Gene name       Gene    Protein_Class
```
### and renamed with 
```bash
{ 
  cat new_header.txt
  tail -n +2 genes_with_protein_type.tsv
} > genes_with_protein_type_cols_renamed.tsv
```
### Then Merged 

this data to my main dataframe with the merge dcript in section 1
```bash
python merge_tsv_by_keys.py \
  --left ranked_genes_with_cluster_categories.tsv \
  --right genes_with_protein_type_cols_renamed.tsv \
  --left-keys "Gene","Gene name" \
  --right-keys "Gene","Gene name" \
  --right-cols "Protein_Class" \
  --out  ranked_genes_with_group_cluster_and_protein_class_data.tsv
```
And then 

### Replaced empty values

with

replace_empty_values.py
```pr

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(
        description="Replace empty or NaN values in a chosen column with a custom value."
    )

    parser.add_argument("-i", "--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV/TSV file")
    parser.add_argument("-c", "--column", required=True, help="Column name to process")
    parser.add_argument("-r", "--replacement", required=True, help="Value to replace empty/NaN")
    parser.add_argument("-d", "--delimiter", default="\t",
                        help="Field delimiter (default: tab for TSV)")

    args = parser.parse_args()

    # Check file exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    print(f"Loading file: {args.input}")
    df = pd.read_csv(args.input, delimiter=args.delimiter)

    if args.column not in df.columns:
        raise ValueError(
            f"Column '{args.column}' not found.\n"
            f"Available columns:\n{list(df.columns)}"
        )

    print(f"Cleaning column: {args.column}")

    # Replace empty strings or whitespace-only strings with NaN
    df[args.column] = df[args.column].replace(r'^\s*$', np.nan, regex=True)

    # Replace NaN with user-supplied replacement value
    df[args.column] = df[args.column].fillna(args.replacement)

    print(f"Saving cleaned file to: {args.output}")
    df.to_csv(args.output, sep=args.delimiter, index=False)

    print("Done!")

if __name__ == "__main__":
    main()
```
### Ran it like
```
python replace_empty_values.py \
    -i ranked_genes_with_group_cluster_and_protein_class_data.tsv \
    -o ranked_genes_cleaned.tsv \
    -c Protein_Class \
    -r Unknown_Protein_Class \
    -d $'\t'
```

## 8-IV Select data

#* Finally I select the number of rows I need in the plot like this

```bash
head -n 10000 final_data_cleaned.tsv> top_10k_from_final_data.tsv
```

# 9) Making interactive plots 

I made interactive plots with universal_plot_maker_plus.py

## 9-I Script

universal_plot_maker_plus.py

```url
https://github.com/TharinduTS/universal_plot_maker_plus/blob/main/README.md
```
## 9-II Run

With following command

```bash
python universal_plot_maker_plus.py \
  --file top_10k_from_final_data.tsv \
  --out Celltype_Enrichment_V2_1_top_10k.html \
  --plot-type bar \
  --x-choices "Gene name | Gene" \
  --y-choices "Enrichment score|log2_enrichment| specificity_tau | Enrichment score (tau penalized)|log2_enrichment_penalized" \
  --default-x "Gene name" \
  --default-y "log2_enrichment_penalized" \
  --color-col "Cell type" \
  --color-choices "Cell type|Cell type group|Cell type class" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class" \
  --title "Celltype Enrichmnt V 2.1" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en
```
# 10) User interface

## 10-I Snapshot

Following is a snapshot of what CellType Enrichment V 2.1 looks like.

<img width="974" height="489" alt="image" src="https://github.com/user-attachments/assets/e8189671-8478-4e14-b6eb-24c61aa00d47" />
A view of user interface

## 10-II Plot behavior 
1	Plot Type – This allows you to select the plot time you want to visualize data in. At the moment, I am only using bar plots as they make most sense for this type of data

2	Color by – This allows you to select which column you want to color your data by. Eg: Cell type gives different colors to different cell types/ Cell type group color bars by the cell type group

3	X – This allows you to select what you want in the X axis. I only gave the options of Gene name and Gene

4	Y- This determines what you see on Y axis

5	Bars to show – This determines how zoomed in/ how many bars you can see in the window. This is specially helpful in exporting a certain number of rows with ‘Export TSV’ button

6	Duplicate policy- This allows you to select how the plot deals with duplicates (Eg: when same Gene is expressed in more than one cell type, how you want them to be visualized. I have explained some use cases below.

7	Primary sort – This lets you select which order you want to sort data by. You have many options with different columns like overall rank by cell type, tau specificity score, gene name etc.

8	Secondary sort- This is used as a tie breaker. When primary sort has the same value for two or more rows, you can use this to decide the sort order.

9	Asc and desc menus- For each sort menus, you can select whether they should be sorted in ascending or descending order. If the column is numeric, these sorts based on value. If the column is string, then it is A-Z or Z-A



## 10-III Filter Menus

1	Cell type class- Lets you filter data by cell type class

2	Cell type group - Lets you filter data by cell type group

3	Cell type - Lets you filter data by cell type

4	Cluster limit – Lets you filter data based on how many clusters were used for the calculation. For now it lets you select all, one cluster or below and 2 clusters or above.

5   Protein_Class - Lets you filter based on protein class

## 10-IV Search Bars

1	Search in Gene – Lets you search genes by gene ID 

2	Search in Gene name – Lets you search genes by gene name

## 10-V Buttons

1	Reset – Resets the plot to original values

2	Export TSV – Exports currently selected data as a TSV file

## 10- VI Information shown

1	Gene – Gene ID

2	Gene name – General name of the gene

3	Cell type - Cell type of currently selected data point/ bar

4	Cell type group - Cell type group of currently selected data point/ bar 

5	Cell type class - Cell type class of currently selected data point/ bar

6	clusters_used – Number of clusters used for currently selected bar/ gene* cell type combination

7	Enrichment score – raw enrichment score

8	log2_enrichment - log2_enrichment score. Good for centering around 0

9	specificity_tau – Specificity score observed by Yanai’s T

10	log2_enrichment_penalized - log2_enrichment score penalized by Yanai’s T

11	top_percent_Cell_type_count – Number of cell types currently selected gene was found in- in the data selected in sections 4 and 6 ( in values more than 0, so it is enrichment/not depletion. And I did not consider lowest 5% of enrichment values in this as I explained In sections 4 and 5)

12	top_percent_Cell_type_group_count - Number of cell type groups currently selected gene was found in- in the data selected in sections 4 and 6 (just like before)

13	top_percent_Cell_type_class_count - Number of cell type classes currently selected gene was found in- in the data selected in sections 4 and 6

14	overall_rank_by_Cell_type – This script ranks gene* cell type combinations (represented by bars in default setting) based on 2 measures, 

•	Number of cell types (or groups or classes) a particular gene is found in (in selected data as explained in section 4 and 6 like before).

•	Enrichment score – In genes that were expressed in a similar number of cell types (or groups or classes, based on sort/filter selection), these combinations are then ranked by enrichment score. 
So basically, if a gene is  only expressed in one cell type, gene with a higher enrichment score is ranked 1 and the gene with lower enrichment score is ranked 2. Any genes that were expressed in more than one cell type comes after both of these genes. 

15	overall_rank_by_Cell_type_group – rank the current selection got just like in overall_rank_by_Cell_type, but this time based on presence in different cell groups instead of different cell types

16	overall_rank_by_Cell_type_class - rank the current selection got just like in overall_rank_by_Cell_type, but this time based on presence in different cell classes instead of different cell types

17	rank_within_Cell_type – rank for current gene based on enrichment score within current celltype

18	rank_within_Cell_type_group - rank for current gene within current cell type group. This does not punish a gene rank for being present in 2 cell types, if they belong to the same cell type group

19	rank_within_Cell_type_class - rank for current gene within current cell type class. This does not punish a gene rank for being present in 2 cell types, if they belong to the same cell type class

20	top_percent_Cell_types – cell types the current gene/s is expressed in (again in the selected data as explained before)

21	top_percent_Cell_type_groups - cell type group/s the current gene is expressed in

22	top_percent_Cell_type_classes - cell type class/es the current gene is expressed in

## 10-VII Use cases

### Visualizing

By Default, CellType Enrichment V2.1 loads with Gene name as X axis and Log2 Enrichment value penalized by specificity score as Y axis. I think this is a good starting point as log2_enrichment puts background level expression level on 0 and it is easier to interpret. 
With defaults set to overlay plots and Primary sort by overall rank by cell type. This plot puts genes that are only expressed in one cell type with highest enrichment score first. Then it plots gradually decreasing enrichment values: still for genes that are only expressed in one cell type (see figure below).

<img width="974" height="300" alt="image" src="https://github.com/user-attachments/assets/a4fbd386-018f-44f9-ab23-8c03b36fff71" />
Overall gene ranking


If you keep scrolling right, then it starts plotting genes that are expressed in two cell types, but highly biased to one cell type. And then genes expressed in two cell types, with a lesser expression bias to one cell type and so on. Different colors represent different cell types here. Figure 4 below shows transition from genes that express in only one cell type to two cell types.

<img width="974" height="370" alt="image" src="https://github.com/user-attachments/assets/c0f063bb-f2b6-491a-b964-997992ea985e" />
Transition from genes expressed in one cell type to genes expressed in two cell types

### Filtering

Filtering Menus gives you the values sorted according to sort columns just like before, but for individual cell type class, cell type group and cell type (see following figure as an example)

<img width="974" height="370" alt="image" src="https://github.com/user-attachments/assets/d235d090-079e-45be-8fc5-5ae16bf8c7c8" />
Transition from genes expressed in one cell type to genes expressed in two cell types - filtered for cell type class blood and immune cells

### Downloading

This script gives you multiple ways to select and download data with ‘Export TSV’ button.

#### Normal download

On the normal view, you can just select and download data for any selection with the tools baked into the plot viewer.

#### A specific number of data points.

Entering a value in the box ‘Bars to show’ lets you select a specific number of data points in the current view and download that number of data points.

#### Sorted data points and top values

In addition to these, you can use sort and filter menus to select your download.

As an example, if you want to download the top values for each cell type, you can set the Primary sort to ‘rank_within_celltype’. This brings all the no 1 ranked genes for each cell type to the top. 

Then you can use secondary sort as ‘cell type’ 

#*NOTE: This last part may need a little more work to get it working perfectly*








