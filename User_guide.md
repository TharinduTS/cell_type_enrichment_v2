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


