# 6) Rank genes on cell specific expression

## 6-I Introduction

rank_genes.py ranks genes by enrichment and summarizes how many groups (e.g., Cell type, Cell type group, Cell type class) each gene appears in within a global top subset. You have full control over:

Grouping column (via --presence-col) used for ranking and “rank within group”. Uniqueness column (via --unique-col) used to count distinct values per gene. Selection & sorting metrics (--top-col, --sorting-col) and top‑% subset size. Output column names, filtering, and NA/negative handling.

This lets you switch between cell‑type specific, group‑level, or class‑level views without changing the code—just alter the CLI flags.

## 6-II CLI help

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

## 6-III Dropping rows closer to background noise

By comparing the enrichment values of known "Cell type specific (CTS)" genes, 

I figured to match the "expression identification sensitivity" of other studies, 

I should not consider the rows with log2_enrichment_penalized value less than 0.5 (with expression values less than about 1.5 times the background expression). Therefore, I started by filtering those rows out.

In other words, I have 167 cell types. In enrichment calculation,

•	If the level of expression of a certain cell specific looking gene that expresses in a certain cell type is ==> X

•	For the rest of the genes, this becomes background expression. If we consider the number of total cell types as 160 for the ease of calculation, background expression level by the highest expressing gene for that cell type becomes ==> X/167

•	The minimum level of expression another gene should have (y) to appear significant to show up for that cell type (to be 1.5 times higher than the base level as I mentioned before) then becomes ==> (X * 1.5)/167 = Y  ==> Y =0.00898X ==> X=100Y(Approximately/not considering the expression of other genes)

•	So in simpler terms, if a certain gene is not expressed at least in a level close to 1% of the level of highly expressing genes for that cell type, that is not considered for the calculation to match the ‘sensitivity level’ of literature I found "Cell type specific (CTS)" genes from. 

•	Then I also have other, drop-na, drop-negatives and top percentage filters I can use, even though it will not make a difference after the filtration step above.

## Filter

Keeps the header, and only rows where log2_enrichment_penalized >= 0.5

```
awk -F'\t' 'NR==1 {print; next} { if ($11+0 >= 0.5) print }' enrich_values_with_cell_class_data.tsv > enrich_values_with_cell_class_data_filtered.tsv
```

## 6-IV Run command

I ran it like following

I Ranked these first by 'Cell type'

#And checked the prescense in the positive enrichment values as minus means depletion. This and drop-na functions should have been already taken care by the previous filtering step. But I am just using them in the first step

I use top-percent value of 100% here as I have performed some hard filtering already

```bash
python rank_genes.py \
  --input enrich_values_with_cell_class_data_filtered.tsv \
  --output ranked_genes_unique_celltypes.tsv \
  --presence-col "Cell type" \
  --unique-col "Cell type" \
  --unique \
  --drop-na \
  --drop-negatives \
  --top-percent 100 \
  --top-col log2_enrichment_penalized \
  --sorting-col log2_enrichment_penalized \
  --gene-col Gene \
  --output-count-col top_percent_Cell_type_count \
  --output-list-col top_percent_Cell_types \
  --output-rank-within-col rank_within_Cell_type \
  --output-overall-rank-col overall_rank_by_Cell_type \
  --verbose
```

Then I rankd them again by the PRESCENCE IN 'Cell type group'

This does no change the previous ranking. This just adds another column to it with ranks based on cell group instead cell type

STILL I DID NOT CALCULATE ENRICHMENT VALUES SEPERATELY FOR 'Cell type group'. ENRICHMENT IS STILL BASED ON ALL THE 'Cell Type' values

AND I AM SETTING --top-percent TO 100% HERE AS IT WAS ALREADY FILTERED IN THE INPUT FILE USED HERE. Same with --drop-na and negatives etc

```bash
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
Then I rankd them again by the PRESCENCE IN 'Cell type class'

```bash
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
