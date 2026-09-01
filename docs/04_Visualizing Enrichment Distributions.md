# 4) Visualizing Enrichment Distributions

## 4-I Introduction

This section provides a simple CLI tool to visualize the distribution of values in any numeric column (e.g., log2_enrichment_penalized) from a TSV/CSV file. It adds:

A horizontal line at 0 Grey shading for the depletion (negative) zone Clear labels for Enrichment and Depletion The ability to highlight near-zero rows with red vertical markers A function to report the % of positive values Optional export of the near-zero rows to a CSV

Why this matters

log₂ > 0 → Enriched (observed higher than expected) log₂ < 0 → Depleted (observed lower than expected) log₂ = 0 → Neutral (observed ≈ expected; enrichment score = 1)

This script makes it easy to see that boundary and inspect data near the neutral point.

## 4-II CLI Help

```txt

usage: plot_distribution.py [-h] -i INPUT -c COLUMN [-d DELIMITER] [-o OUTPUT]
                            [--near-zero NEAR_ZERO] [--export-near-zero EXPORT_NEAR_ZERO]

Visualize distribution and analyze enrichment values.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input CSV/TSV file
  -c COLUMN, --column COLUMN
                        Column name to plot distribution
  -d DELIMITER, --delimiter DELIMITER
                        Delimiter (default: tab)
  -o OUTPUT, --output OUTPUT
                        Output plot filename
  --near-zero NEAR_ZERO
                        Show N rows closest to zero
  --export-near-zero EXPORT_NEAR_ZERO
                        Export near-zero rows to CSV file
```

4-III Run

I ran it like following

```bash
python plot_distribution.py -i enrichment_values_for_filtered_celltypes.tsv -c log2_enrichment_penalized --near-zero 10 -o log2_enrichment_penalized_distribution.png
```

This code plots the distribution, Give you the top percentage of rows with positive values and shows you the rows around that value

```bash
Percentage of rows with positive log2_enrichment_penalized: 18.17%
```

With this data, in the ranking step below (section 5), I am

	1.dropping all the negative values and 0s as they represent cell type*gene combinations with equal or less than expression to the background
	2.Checking the presense of expression of genes in cell types in top x% of the enrichment values FROM THE REMAINING values so I only get gene expression which is significantly higher than background
log2_enrichment_penalized_distribution.png

<img width="1000" height="700" alt="log2_enrichment_penalized_distribution" src="https://github.com/user-attachments/assets/ad770994-6ef6-4b93-a0b4-f0fe4915e96b" />
