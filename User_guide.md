# Cell type enrichment User Guide

## User interface

=======
Following is a snapshot of Cell Type Enrichment V.2 interface

<img width="1326" height="941" alt="User_interface_screenshot" src="https://github.com/user-attachments/assets/670f8e6d-a1b9-4e88-aedb-4af3eee8f469" />

## Menus and Filters
### 1) Plot behavior
1 Plot Type – This allows you to select the plot time you want to visualize data in. At the moment, I am only using bar plots as they make most sense for this type of data

2 Color by – This allows you to select which column you want to color your data by. Eg: Cell type gives different colors to different cell types/ Cell type group color bars by the cell type group

3 X – This allows you to select what you want in the X axis. I only gave the options of Gene name and Gene

4 Y- This determines what you see on Y axis

5 Bars to show – This determines how zoomed in/ how many bars you can see in the window. This is specially helpful in exporting a certain number of rows with ‘Export TSV’ button

6 Duplicate policy- This allows you to select how the plot deals with duplicates (Eg: when same Gene is expressed in more than one cell type, how you want them to be visualized. I have explained some use cases below.

7 Primary sort – This lets you select which order you want to sort data by. You have many options with different columns like overall rank by cell type, tau specificity score, gene name etc.

8 Secondary sort- This is used as a tie breaker. When primary sort has the same value for two or more rows, you can use this to decide the sort order.

9 Asc and desc menus- For each sort menus, you can select whether they should be sorted in ascending or descending order. If the column is numeric, these sorts based on value. If the column is string, then it is A-Z or Z-A

### 2) Filter Menus
1 Cell type class- Lets you filter data by cell type class

2 Cell type group - Lets you filter data by cell type group

3 Cell type - Lets you filter data by cell type

4 Cluster limit – Lets you filter data based on how many clusters were used for the calculation. Currently, it lets you select all, one cluster or below and 2 clusters or above.

### 3) Search Bars
1 Search in Gene – Lets you search genes by gene ID

2 Search in Gene name – Lets you search genes by gene name

3 Search in Present tissues - Lets you search gene / cell type combinations by expressed tissue type

### 4) Buttons
1 Reset – Resets the plot to original values

2 Export TSV – Exports currently selected data as a TSV file

### 5) Information shown
1 Gene – Gene ID

2 Gene name – General name of the gene

3 Cell type - Cell type of currently selected data point/ bar

4 Cell type group - Cell type group of currently selected data point/ bar

5 Cell type class - Cell type class of currently selected data point/ bar

6 clusters_used – Number of clusters used for currently selected bar/ gene* cell type combination

7 Enrichment score – raw enrichment score

8 log2_enrichment - log2_enrichment score. Good for centering around 0

9 specificity_tau – Specificity score observed by Yanai’s T

10 log2_enrichment_penalized - log2_enrichment score penalized by Yanai’s T

11 top_percent_Cell_type_count – Number of cell types currently selected gene was found in- in the data selected in sections 4 and 6 ( in values more than 0, so it is enrichment/not depletion. And I did not consider lowest 5% of enrichment values in this as I explained In sections 4 and 5)

12 top_percent_Cell_type_group_count - Number of cell type groups currently selected gene was found in- in the data selected in sections 4 and 6 (just like before)

13 top_percent_Cell_type_class_count - Number of cell type classes currently selected gene was found in- in the data selected in sections 4 and 6

14 overall_rank_by_Cell_type – This script ranks gene* cell type combinations (represented by bars in default setting) based on 2 measures,

• Number of cell types (or groups or classes) a particular gene is found in (in selected data as explained in section 4 and 6 like before).

• Enrichment score – In genes that were expressed in a similar number of cell types (or groups or classes, based on sort/filter selection), these combinations are then ranked by enrichment score. So basically, if a gene is only expressed in one cell type, gene with a higher enrichment score is ranked 1 and the gene with lower enrichment score is ranked 2. Any genes that were expressed in more than one cell type comes after both of these genes.

15 overall_rank_by_Cell_type_group – rank the current selection got just like in overall_rank_by_Cell_type, but this time based on presence in different cell groups instead of different cell types

16 overall_rank_by_Cell_type_class - rank the current selection got just like in overall_rank_by_Cell_type, but this time based on presence in different cell classes instead of different cell types

17 rank_within_Cell_type – rank for current gene based on enrichment score within current celltype

18 rank_within_Cell_type_group - rank for current gene within current cell type group. This does not punish a gene rank for being present in 2 cell types, if they belong to the same cell type group

19 rank_within_Cell_type_class - rank for current gene within current cell type class. This does not punish a gene rank for being present in 2 cell types, if they belong to the same cell type class

20 top_percent_Cell_types – cell types the current gene/s is expressed in (again in the selected data as explained before)

21 top_percent_Cell_type_groups - cell type group/s the current gene is expressed in

22 top_percent_Cell_type_classes - cell type class/es the current gene is expressed in

23 Present tissues - Tissue(s) current gene / cell type combination shows expression in



