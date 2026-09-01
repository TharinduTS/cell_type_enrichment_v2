
#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_positive_percentage(series):
    total = len(series)
    positive_count = (series > 0).sum()
    return (positive_count / total) * 100

def get_rows_near_zero(df, column, n):
    # Drop NA values and sort by absolute distance from zero
    df_clean = df.loc[df[column].notna()].copy()
    df_clean["abs_val"] = df_clean[column].abs()
    return df_clean.sort_values("abs_val").drop(columns="abs_val").head(n)

def main():
    parser = argparse.ArgumentParser(description="Visualize distribution and analyze enrichment values.")
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV/TSV file")
    parser.add_argument("-c", "--column", required=True, help="Column name to plot distribution")
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter (default: tab)")
    parser.add_argument("-o", "--output", default=None, help="Output plot filename")
    parser.add_argument("--near-zero", type=int, default=0, help="Show N rows closest to zero")
    parser.add_argument("--export-near-zero", default=None, help="Export near-zero rows to CSV file")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input, delimiter=args.delimiter)

    if args.column not in df.columns:
        print(f"Error: Column '{args.column}' not found. Available columns: {list(df.columns)}")
        return

    # Positive percentage
    pct_positive = calculate_positive_percentage(df[args.column])
    print(f"Percentage of rows with positive {args.column}: {pct_positive:.2f}%")

    # Plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 7))
    ax = sns.histplot(df[args.column], kde=True, color="blue", bins=30)

    # Vertical line at 0
    plt.axvline(0, color="black", linestyle="--", linewidth=1)

    # Get axis limits after plotting
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    # Shade depletion zone
    plt.axvspan(xmin, 0, color="grey", alpha=0.2)

    # Labels for zones
    plt.text(xmin + (xmax - xmin) * 0.05, ymax * 0.92, "Depletion", color="grey", fontsize=12)
    plt.text(xmin + (xmax - xmin) * 0.70, ymax * 0.92, "Enrichment", color="blue", fontsize=12)

    # Highlight near-zero rows if requested
    if args.near_zero > 0:
        near_zero_rows = get_rows_near_zero(df, args.column, args.near_zero)
        for val in near_zero_rows[args.column]:
            plt.axvline(val, color="red", linestyle=":", linewidth=1)
        print(f"\nTop {args.near_zero} rows closest to zero in {args.column}:")
        print(near_zero_rows)

        # Export if requested
        if args.export_near_zero:
            near_zero_rows.to_csv(args.export_near_zero, index=False)
            print(f"Near-zero rows exported to {args.export_near_zero}")

    # Title and labels
    plt.title(f"Distribution of {args.column}")
    plt.xlabel(args.column)
    plt.ylabel("Frequency")

    # Save plot
    output_file = args.output if args.output else f"{args.column}_distribution.png"
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved as {output_file}")

if __name__ == "__main__":
    main()