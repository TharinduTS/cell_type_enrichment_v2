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

23 Protein class - Any predicted protein classes

24 Tissues - Tissue(s) current gene / cell type combination shows expression in

## Use cases

### Visualizing
By Default, CellType Enrichment V2 loads with Gene name as X axis and Log2 Enrichment value penalized by specificity score as Y axis. I think this is a good starting point as log2_enrichment puts background level expression level on 0 and it is easier to interpret. With defaults set to overlay plots and Primary sort by overall rank by cell type. This plot puts genes that are only expressed in one cell type with highest enrichment score first. Then it plots gradually decreasing enrichment values: still for genes that are only expressed in one cell type (see figure below).

<img width="1337" height="427" alt="initial_plot_image" src="https://github.com/user-attachments/assets/0288e100-ee2c-4690-b2f5-5cf3ae3dac3d" />

#### Tissue specific expression subplot

Once you click a bar representing a gene, this will automatically open a subplot on the side, showing how the expression of that gene changes in different tissue types.

<img width="365" height="312" alt="tissue plot" src="https://github.com/user-attachments/assets/4375e244-942e-4fb7-aa39-1642ab3cfa45" />


### Overall gene ranking
If you keep scrolling right, then it starts plotting genes that are expressed in two cell types, but highly biased to one cell type. And then genes expressed in two cell types, with a lesser expression bias to one cell type and so on. Different colors represent different cell types here. Figure 4 below shows transition from genes that express in only one cell type to two cell types.

<img width="1281" height="365" alt="scrolled_look" src="https://github.com/user-attachments/assets/889c10bf-5e22-45ed-8c0b-95519428128f" />

### Filtering
Filtering Menus gives you the values sorted according to sort columns just like before, but for individual cell type class, cell type group and cell type

### Downloading
#NOTE #When you download rows based on a order number, make sure you have selected separate as duplicate policy so the overlay/stacked plots for the same gene do not confuse you

This script gives you multiple ways to select and download data with ‘Export TSV’ button.

#### Normal download
On the normal view, you can just select and download data for any selection with the tools baked into the plot viewer.

#### A specific number of data points.
Entering a value in the box ‘Bars to show’ lets you select a specific number of data points in the current view and download that number of data points.

#### Sorted data points and top values
In addition to these, you can use sort and filter menus to select your download.

As an example, if you want to download the top values for each cell type, you can set the Primary sort to ‘rank_within_celltype’. This brings all the no 1 ranked genes for each cell type to the top.

Then you can use secondary sort as ‘cell type’

## App Walkthrough
In addition to all the features I had explained in my previous version, following are a few of the new, interesting features I thought worth noting. In this section, I try to explain these with examples.

### Ranking by cell type groups and cell type classes
In addition to the default loading setting: genes ranked by cell type specificity, now this supports genes ranked by cell type groups and classes. As an example, most germ cells do not make it to the top of the list when you are in the default ‘overall rank by cell type’ setting as they are expressed in multiple cell types that are closely related. But if you change it to ‘overall rank by cell type group’, this does not push a gene back for being present in multiple cell types- if they are closely related (See figure : Gene SPEM2).

<img width="975" height="370" alt="image" src="https://github.com/user-attachments/assets/bb3ffc3a-5d6b-4f65-8305-d2e702bdff81" />
Overall rank by cell type group brings germ cells to the top of the list

This gene SPEM2 in the example is present in early spermatids and late spermatids, which sometimes you would not want to consider as two different cell types. This brings such genes up in the list. You can do the same with cell type group (please note that It is still using enrichment values calculated by cell type for now. This works totally fine even if this is not the 100% perfect solution). 

### Color by groups

Sometimes the amount of details can be too much depending on what you are looking for. Therefore, now you have properly working color by options.
Using the same example above, if you change the color by option to ‘Cell type group’, now you can see the whole bar for the gene SPEM2 and genes in the same group (like PRM2 right next to SPEM2 in figure 6) have the same color.

<img width="975" height="459" alt="image" src="https://github.com/user-attachments/assets/20df0fa3-15db-46a5-88b9-ee4f3ef2a6b9" />
Color by cell type group bringing cell types from the same group together

Again, you can do this by class as well.

### Duplicate policies

After the top half of the list, genes can be expressed in more than one cell type. Duplicate policies give you multiple ways of handling this.

#### Overlay

The default setting: overlay plots the cell types that express the same gene together. This is like a 3D plot with the cell type with the highest enrichment value in the back and the cell types that express the same gene, but with lower enrichment values in the front. You can see this in the same example (please see figure 8). SPEM2 expression in late spermatids (represented by light pink) is higher than early spermatids (dark pink bar in front of the previous).

<img width="338" height="455" alt="image" src="https://github.com/user-attachments/assets/a033cee7-69b0-49e8-986f-c3bb9d45c6fc" />

Overlay plots showing how the same gene is expressed in multiple cell types

#### Stack

This option lets you stack enrichment values for all the cell types for a certain gene on top of each other. This is good for broad level comparisons.

#### Sum

This gives you the sum of enrichment values for each gene for all the cell types it is expressed in. I found this nice when you compare expression of different cell type groups, and you are not particularly interested in some sub cell types.  Using the same example above, even though it looked like the top ranking gene ‘LALBA’ had a way higher level of cell specific expression in the default setting, when you check the sum of the expression on SPEM2 in closely related early and late spermatids, the difference noticeably decreases 

<img width="975" height="471" alt="image" src="https://github.com/user-attachments/assets/496c02ec-2acc-4408-b8ba-2a043a756b64" />
Using duplicate policy 'sum' for cell group comparisons

#### Separate

Even though previous policies are nice for visualizing, they can be confusing when you try to download data based on their rank. Using the same example, the cell types SPEM2 is expressed in, early and late spermatids are plotted together in previous modes despite their overall rank by cell type groups being 245 and 452. Separate mode splits such groupings and adds numbers to genes that are expressed in more than one cell type to show the rank for that particular cell type within that gene. See the example in figures 10 and 11 showing how the ranking is visualized for GGTLC1 expressed in alveolar cells type 1 (Red) and alveolar cells type 2 (Purple) in ‘overlay’ and ‘separate’ plots.

<img width="960" height="520" alt="image" src="https://github.com/user-attachments/assets/d1f6db8c-843e-42a7-81d9-5c6d26576d3e" />

Genes expressed in two cell types represented in 'overlay' plots

<img width="975" height="574" alt="image" src="https://github.com/user-attachments/assets/70c8ae0e-846d-4f57-ac63-21a57f0387a5" />

Genes expressed in two cell types represented in 'separate' plots

#### Max, Mean, Median, First

Gives you option to see Max, Mean, Median or first occurrence

### Protein Classes

In addition to previous filtering tools, I added a protein class filter with following classes

Predicted intracellular proteins

Predicted membrane proteins

Predicted secreted proteins

According to HPA (https://www.proteinatlas.org/humanproteome/tissue/secretome#classification_of_the_secretome), The human secretome comprises all proteins that are potentially secreted from the cell.






