## 8-II Merge protein type data

I Started by renaming the columns of genes_with_protein_type.tsv in a meaningful way to make it easier to map

First I

made a mapping tsv
file for column names with new column names

new_header.txt

```txt
Gene name       Gene    Protein_Class
```

#MAKE SURE THESE HEADERS ARE TAB SEPERATED. SOMETIMES COPY PASTING DOES NOT WORK

and renamed with

```bash
{ 
  cat new_header.txt
  tail -n +2 genes_with_protein_type.tsv
} > genes_with_protein_type_cols_renamed.tsv
```

Then Merged
this data to my main dataframe with the merge dcript in section 1

```
python merge_tsv_by_keys.py \
  --left ranked_genes_with_cluster_categories.tsv \
  --right genes_with_protein_type_cols_renamed.tsv \
  --left-keys "Gene","Gene name" \
  --right-keys "Gene","Gene name" \
  --right-cols "Protein_Class" \
  --out  ranked_genes_with_group_cluster_and_protein_class_data.tsv
```

And then

Replaced empty values
with

replace_empty_values.py

### Ran it like

```py
python replace_empty_values.py \
    -i ranked_genes_with_group_cluster_and_protein_class_data.tsv \
    -o ranked_genes_cleaned.tsv \
    -c Protein_Class \
    -r No_predicted_protein_class \
    -d $'\t'
```