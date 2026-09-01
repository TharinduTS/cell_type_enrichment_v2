Merge Script CLI help

#************** NOTE ************************************

This script can leave rows without matching keys on right table column value empty. Therefore you should use a command like following to find the rows with unassigned values to avoid errors down the pipeline (I have dealt with this issue in chapter 5) see the example below from section 5

```
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Cell type group") col=i; next} $col=="" {print}' enrich_values_with_cell_class_data.tsv | cut -f 3 | sort -u
```


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
