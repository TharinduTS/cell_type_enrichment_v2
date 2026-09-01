# 5) Add celltype info layers

## 5-I Introduction

My enrichment value table only has cell type information. But it can be useful to have information like cell type group and cell type class. I am adding those layers to the dataset here

HPA has these information on the table here:

```http
https://www.proteinatlas.org/download/tsv/rna_single_cell_type_cell_types.tsv.zip
```
I can use the same merging script from section 1

## 5-II Run command

Run it like following

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

```txt
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
