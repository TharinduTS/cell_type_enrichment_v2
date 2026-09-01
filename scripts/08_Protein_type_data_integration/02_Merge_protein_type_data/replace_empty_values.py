#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(
        description="Replace empty or NaN values in a chosen column with a custom value."
    )

    parser.add_argument("-i", "--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV/TSV file")
    parser.add_argument("-c", "--column", required=True, help="Column name to process")
    parser.add_argument("-r", "--replacement", required=True, help="Value to replace empty/NaN")
    parser.add_argument("-d", "--delimiter", default="\t",
                        help="Field delimiter (default: tab for TSV)")

    args = parser.parse_args()

    # Check file exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    print(f"Loading file: {args.input}")
    df = pd.read_csv(args.input, delimiter=args.delimiter)

    if args.column not in df.columns:
        raise ValueError(
            f"Column '{args.column}' not found.\n"
            f"Available columns:\n{list(df.columns)}"
        )

    print(f"Cleaning column: {args.column}")

    # Replace empty strings or whitespace-only strings with NaN
    df[args.column] = df[args.column].replace(r'^\s*$', np.nan, regex=True)

    # Replace NaN with user-supplied replacement value
    df[args.column] = df[args.column].fillna(args.replacement)

    print(f"Saving cleaned file to: {args.output}")
    df.to_csv(args.output, sep=args.delimiter, index=False)

    print("Done!")

if __name__ == "__main__":
    main()