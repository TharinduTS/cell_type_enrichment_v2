# 9) Adding tissue expression profile

Introduction:
In addition to having all the information layers above, it would be nice to see how these genes with cell type biased expression are expressed in different tissue types. Therefore I am trying to extract and add this data layer and visualize that in the same page, but without disturbing the original clean layout.

9-I Tissue expression layer
I am going to use the same enrichment calculation script I used in chapter 3, but this time with tissue data layer.

I started by copying my filtered data file and enrichment calculation script to a new directory.

First I fused Tissue and Cell type columns
because I need both those data layers considered when I decide categories for enrichment calculation

```bash
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print $0,"Tissue:CellType"} NR>1{print $0,$1":"$5}' combined_expression_data_filtered.tsv > combined_expression_data_with_tissue_celltype_fused.tsv
```

calculated enrichment values
with this newly fused column for tissue:cell type combinations

```bash
./run_celltype_enrichment_v1_4.sh --input-file combined_expression_data_with_tissue_celltype_fused.tsv --output-file enrichment_values_for_tissue_types.tsv --min-clusters 1 --min-count 50 --specificity-mode penalize --min-specificity 1 --cell-type-col "Tissue:CellType"  --batch-col "Tissue:CellType"
```

#Note how I added two new commands here

```bash
--cell-type-col "Tissue:CellType"  --batch-col "Tissue:CellType"
```

defining which column should be used for grouping

This gives me enrichment values for Tissue:cellType combinations like following

Note that the progress messages given by this script can be confusing now as this script was written calculating enrichment values for cell types in mind.Cell types now mean Tissue:CellType combinations

```txt
Gene    Gene name       Tissue:CellType avg_nCPM        weight_sum      clusters_used   specificity_tau Enrichment score        log2_enrichment Enrichment score (tau penalized)        log2_enrichment_penalized       single_cell_type_gene
ENSG00000167531 LALBA   breast:breast lactating cells   483428.6        420     1       0.9999997849673347      4606146.319293296       22.135128809974248      4606145.328821376       22.135128499747655      False
ENSG00000135222 CSN2    breast:breast lactating cells   333043.0        420     1       0.9999991703984089      1201050.8623900507      20.195865817257094      1201049.8659963442      20.1958646203945        False
```

selected only the useful column for clarity

```bash
less enrichment_values_for_tissue_types.tsv | cut -f 1,2,3,6,11 > selected_tissue_expression_columns.tsv
```
### splitted this combined column
into 2 columns again for ease of data handling

```awk
awk -F'\t' '
BEGIN { OFS="\t" }

NR==1 {
    # Print new header
    print $1, $2, "Tissue", "Cell type", $4, $5
    next
}

{
    split($3, a, ":")
    print $1, $2, a[1], a[2], $4, $5
}
' selected_tissue_expression_columns.tsv> combined_expression_data_split.tsv
```
Merged cluster info into Tissue data

```awk
awk -F'\t' '
BEGIN { OFS="\t" }
NR==1 { print; next }
{
    $3 = $3 " with " $5
    print
}
' combined_expression_data_split.tsv > combined_expression_data_split_with_clusters.tsv
```
## 9-II Merge data
Then I added this data to my main dataframe. To avoid further complicating the way this data is shown, I merged this data matching columns "Gene Gene name Cell type" and looking for the Tissue inside the 'Present Tissues' column with following script

Script
add_enrichment_to_tissues.py

### CLI help

```txt
add_enrichment_to_tissues.py
============================

Propagate enrichment values from one TSV table into multi-tissue annotations
in another TSV file.

For each row in Table 2, the script:
- Matches rows from Table 1 using: Gene, Gene name, and Cell type
- Checks whether the Tissue from Table 1 appears inside the "Present tissues" column
- Prepends the enrichment value to each tissue as: value:tissue
- Assigns a user-defined value to tissues with no matching enrichment

USAGE
-----
python3 add_enrichment_to_tissues.py \
  --table1 TABLE1.tsv \
  --table2 TABLE2.tsv \
  --output OUTPUT.tsv \
  [options]

REQUIRED ARGUMENTS
------------------
--table1     TSV file containing tissue-specific enrichment values
--table2     TSV file containing multi-tissue annotations
--output     Output TSV file

OPTIONAL ARGUMENTS
------------------
--gene-col                 Gene ID column (default: Gene)
--gene-name-col            Gene name column (default: Gene name)
--cell-type-col            Cell type column (default: Cell type)

--tissue-col-table1        Tissue column in table1 (default: Tissue)
--tissue-col-table2        Multi-tissue column in table2 (default: Present tissues)
--enrichment-col           Value to propagate (default: log2_enrichment_penalized)
--tissue1-trim-after "-"   Ignores the part after this separator in matching
--label-source table1      This lets you keep the value after the trim character above in the output/ leaving cluster data in the output, even though not used in match

--tissue-sep               Tissue delimiter (default: &)
--value-sep                Separator between value and tissue (default: :)
--missing-value            Value for tissues without enrichment (default: NA)

--case-insensitive         Enable case-insensitive matching
--output-tissues-col       Write enriched output to a new column instead of overwriting

EXAMPLE
-------
python3 add_enrichment_to_tissues.py \
  --table1 combined_expression_data_split.tsv \
  --table2 ranked_genes_cleaned_with_tissue_data.tsv \
  --output Final_data_with_tissue_expression_data.tsv \
  --missing-value Not_high_enough \
  --tissue-sep "&" \
  --value-sep ":" \
  --case-insensitive

EXAMPLE TRANSFORMATION
----------------------
Input:
  stomach & tongue & breast

Output:
  Not_high_enough:stomach & Not_high_enough:tongue & 22.13:breast

NOTES
-----
- Every tissue always receives a value
- Original row structure is preserved
- Suitable for large TSV gene expression datasets
```

### Duplicate tissue data
Because the next step is going to add extra data to the Tissue column, I am keeping a copy of that column to be used in 'color by tissue' later in plotting

```awk
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print $0, "Tissues"; next} {print $0, $29}' ranked_genes_cleaned_with_tissue_data.tsv > ranked_genes_cleaned_with_tissue_data_duplicated.tsv
```

### Run command

```bash
python3 add_enrichment_to_tissues.py \
  --table1 combined_expression_data_split_with_clusters.tsv \
  --table2 ranked_genes_cleaned_with_tissue_data_duplicated.tsv \
  --output Final_data_with_tissue_expression_data.tsv \
  --missing-value "Not_high_enough" \
  --value-sep ":" \
  --tissue-sep "&" \
  --tissue1-trim-after " with " \
  --label-source table1 \
  --case-insensitive
```

9-III Select data to plot
Just like I did in chapter 8, here I am selecting how many top data rows to plot, to avoid procssing complications

```bash
head -n 50000 Final_data_with_tissue_expression_data.tsv> top_50k_from_final_data.tsv
```
11-IV Making interactive plots
Then I made interactive plots with subplots to show tissue expression profile with the improved version of universal_plot_maker to include sub plots.

script
Updated universal_plot_maker_plus.py can be found here

```http
https://github.com/TharinduTS/universal_plot_maker_plus_with_subplot/blob/main/README.md
```


### Run command
I ran it with the selected data from before, with following command

50k data points for managable data chunks

```bash
python universal_plot_maker_plus.py \
  --file top_50k_from_final_data.tsv \
  --out Celltype_Enrichment_V2_2_top_50k.html \
  --plot-type bar \
  --x-choices "Gene name | Gene" \
  --y-choices "Enrichment score|log2_enrichment| specificity_tau | Enrichment score (tau penalized)|log2_enrichment_penalized" \
  --default-x "Gene name" \
  --default-y "log2_enrichment_penalized" \
  --color-col "Cell type" \
  --color-choices "Cell type|Cell type group|Cell type class|Tissues" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name|Present tissues" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class|Tissues" \
  --title "Celltype Enrichmnt V 2.2" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en \
  --pt-enable \
  --pt-col "Present tissues" \
  --pt-title "Enrichment per present tissue" \
  --pt-x-label "Tissue" \
  --pt-y-label "log2 Enrichment Penalized" \
  --pt-color "#2a9d8f" \
  --pt-height 360 \
  --pt-width auto \
  --pt-rotate -35 \
  --pt-container-id "present-tissues-plot" \
  --pt-enable --pt-mode flow \
  --pt-anchor "#rowDetails" --pt-position after \
  --pt-offset-x -300 --pt-offset-y -10
```
10k data points for emailing

```bash
python universal_plot_maker_plus.py \
  --file top_10k_from_final_data.tsv \
  --out Celltype_Enrichment_V2_2_top_10k.html \
  --plot-type bar \
  --x-choices "Gene name | Gene" \
  --y-choices "Enrichment score|log2_enrichment| specificity_tau | Enrichment score (tau penalized)|log2_enrichment_penalized" \
  --default-x "Gene name" \
  --default-y "log2_enrichment_penalized" \
  --color-col "Cell type" \
  --color-choices "Cell type|Cell type group|Cell type class|Tissues" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name|Present tissues" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class|Tissues" \
  --title "Celltype Enrichmnt V 2.2" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en \
  --pt-enable \
  --pt-col "Present tissues" \
  --pt-title "Enrichment per present tissue" \
  --pt-x-label "Tissue" \
  --pt-y-label "log2 Enrichment Penalized" \
  --pt-color "#2a9d8f" \
  --pt-height 360 \
  --pt-width auto \
  --pt-rotate -35 \
  --pt-container-id "present-tissues-plot" \
  --pt-enable --pt-mode flow \
  --pt-anchor "#rowDetails" --pt-position after \
  --pt-offset-x -300 --pt-offset-y -10
```




