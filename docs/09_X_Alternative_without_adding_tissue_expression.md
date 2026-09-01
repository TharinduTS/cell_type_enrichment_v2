# 9) Making interactive plots
I made interactive plots with universal_plot_maker_plus.py

## 9-I Script
universal_plot_maker_plus.py

## 9-II Run
I made the plots with half the data to make it easily sharable (to keep it under 25mb) With following command

### Top 10k rows - For sample emailing

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
  --color-choices "Cell type|Cell type group|Cell type class|Present tissues" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name|Present tissues" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class|Present tissues" \
  --title "Celltype Enrichmnt V 2.2" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en
```
## Top 50k rows - For all important data points
This also takes not putting too much stress on the browser and keeping it easily usable without lag into account


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
  --color-choices "Cell type|Cell type group|Cell type class|Present tissues" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name|Present tissues" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class|Present tissues" \
  --title "Celltype Enrichmnt V 2.2" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en
```
And made plots with all the remaining data after filtering with the following command

All rows

```
python universal_plot_maker_plus.py \
  --file ranked_genes_cleaned_with_tissue_data.tsv \
  --out Celltype_Enrichment_V2_2.html \
  --plot-type bar \
  --x-choices "Gene name | Gene" \
  --y-choices "Enrichment score|log2_enrichment| specificity_tau | Enrichment score (tau penalized)|log2_enrichment_penalized" \
  --default-x "Gene name" \
  --default-y "log2_enrichment_penalized" \
  --color-col "Cell type" \
  --color-choices "Cell type|Cell type group|Cell type class|Present tissues" \
  --filter-cols "Cell type class|Cell type group|Cell type|cluster_limit|Protein_Class" \
  --search-cols "Gene|Gene name|Present tissues" \
  --details "Gene|Gene name|Cell type|Cell type group|Cell type class|clusters_used|Enrichment score|log2_enrichment| specificity_tau |log2_enrichment_penalized|top_percent_Cell_type_count|top_percent_Cell_type_group_count|top_percent_Cell_type_class_count|overall_rank_by_Cell_type|overall_rank_by_Cell_type_group|overall_rank_by_Cell_type_class|rank_within_Cell_type|rank_within_Cell_type_group|rank_within_Cell_type_class|top_percent_Cell_types|top_percent_Cell_type_groups|top_percent_Cell_type_classes|Protein_Class|Present tissues" \
  --title "Celltype Enrichmnt V 2.2" \
  --dup-policy overlay \
  --sort-primary "overall_rank_by_Cell_type" \
  --sort-primary-order asc \
  --sort-secondary "log2_enrichment_penalized" \
  --sort-secondary-order desc \
  --initial-zoom 100 \
  --self-contained \
  --lang en
```
