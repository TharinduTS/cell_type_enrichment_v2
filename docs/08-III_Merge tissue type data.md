## 8-III Merge tissue type data

Then I added tissue data for these filtered rows from the earlier filtered data file "combined_expression_data_filtered.tsv" from chapter 2

Script
match_and_collect.py

### CLI help
Match_and_collect.py matches combinations of columns between two tables, collects values from the second table, and appends them as a new column to the first table. Supports CSV, TSV, TXT, XLSX, and XLS files, flexible filtering, deduplication, sorting, and case-insensitive matching.


Usage

```bash
python match_and_collect.py [OPTIONS] file1 file2
```

```txt

match_and_collect.py — CLI Help
================================

Match combinations of columns between two tables, collect values from the second, and append them to the first.
Supports CSV/TSV/TXT/XLSX/XLS, flexible filtering, deduplication, sorting, and case-insensitive matching.

USAGE
-----
python match_and_collect.py [OPTIONS] file1 file2

REQUIRED
--------
file1                         Path to the first file (CSV/TSV/TXT/XLSX/XLS). Receives the new column.
file2                         Path to the second file (CSV/TSV/TXT/XLSX/XLS). Source of collected values.
--file1-keys COL [COL ...]    Key columns in file1 (e.g., "Gene name" "Cell type").
--file2-keys COL [COL ...]    Corresponding key columns in file2 (same count/order as file1).
--collect-column COL          Column from file2 to collect (e.g., "Tissue").
--output-column NAME          Name of the new column added to file1.

FILTERING (choose ONE method)
-----------------------------
Simple:
  --filter-col COL            Column to filter on in file2 (e.g., nCPM).
  --filter-op OP              Operator: > >= < <= == !=
  --filter-value VAL          Value for comparison (e.g., 0)
Advanced (pandas query):
  --filter-expr "EXPR"        e.g., "`nCPM` > 0 and `Cell count` >= 10"
(Do not mix simple and advanced filtering.)

NORMALIZATION
-------------
--case-insensitive            Lowercase key values before matching.
--no-strip                    Do not trim whitespace in key columns (default trims).

AGGREGATION & OUTPUT FORMAT
---------------------------
--sep-out SEP                 Separator for collected values (default: "; ").
--drop-dupes                  Remove duplicate collected values.
--sort-values                 Sort collected values alphabetically.

FILE FORMAT SETTINGS
--------------------
--file1-sep SEP               Override delimiter for file1 (CSV/TXT/TSV).
--file2-sep SEP               Override delimiter for file2 (CSV/TXT/TSV).
--sheet1 NAME                 Excel sheet name for file1 (.xlsx/.xls).
--sheet2 NAME                 Excel sheet name for file2 (.xlsx/.xls).

RESULT WRITING
--------------
--output PATH                 Output file (CSV/TSV/TXT/XLSX/XLS). If omitted, prints a preview only.

HELP
----
-h, --help                    Show this help and exit.

EXAMPLES
--------
1) Basic match on two keys; collect Tissue where nCPM > 0; write Excel:
  python match_and_collect.py table1.xlsx table2.tsv \
    --file1-keys "Gene name" "Cell type" \
    --file2-keys "Gene name" "Cell type" \
    --collect-column "Tissue" \
    --filter-col nCPM --filter-op ">" --filter-value 0 \
    --output-column "Matched tissues" \
    --sep-out "; " --drop-dupes \
    --output merged.xlsx

2) Different key names; advanced filter; write CSV:
  python match_and_collect.py table1.csv table2.csv \
    --file1-keys "Gene name" "Cell type" \
    --file2-keys Gene "Cell type" \
    --collect-column Tissue \
    --filter-expr "`nCPM` > 0 and `Cell count` >= 10" \
    --output-column "Filtered tissues" \
    --output out.csv

NOTES
-----
- Key column lists must match in length and order.
- Quote column names with spaces. In --filter-expr, wrap spaced names in backticks.
- TSV auto-detected by extension; otherwise set --file1-sep and/or --file2-sep.
```

### Run command

```bash
python match_and_collect.py \
  ranked_genes_cleaned.tsv combined_expression_data_filtered.tsv \
  --file1-keys "Gene name" "Cell type" \
  --file2-keys "Gene name" "Cell type" \
  --collect-column "Tissue" \
  --filter-col nCPM --filter-op ">" --filter-value 0 \
  --output-column "Present tissues" \
  --sep-out " & " \
  --drop-dupes \
  --output ranked_genes_cleaned_with_tissue_data.tsv
```

## 8-IV Select data

Finally I select the number of rows I need in the plot like this

```bash
head -n 10000 ranked_genes_cleaned_with_tissue_data.tsv> top_10k_from_final_data.tsv
```
