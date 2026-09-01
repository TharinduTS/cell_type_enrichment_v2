# 7) Cluster data integration

## 7-I Introduction

#*First I added cluster categories (to add as a filter later on section 8) with a minimum limit with this script. You can select a custom number of clusters to filter by. Here I use

1 or below and

2 or above as limits

## 7-II CLI help

```bash
python categorize_column.py <input_file> <output_file> <selected_column> <selected_number> <result_column>
```

7-III Run command

I ran it like following

```bash
python categorize_column.py ranked_genes_unique_celltypes_and_groups_and_classes.tsv ranked_genes_with_cluster_categories.tsv clusters_used 2 cluster_limit
```