import pandas as pd
import sys

def main():
    if len(sys.argv) != 6:
        print("Usage: python categorize_column.py <input_file> <output_file> <selected_column> <selected_number> <result_column>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    selected_column = sys.argv[3]
    selected_number = float(sys.argv[4])
    result_column = sys.argv[5]

    # Load input file (auto-detect delimiter)
    try:
        df = pd.read_csv(input_file, sep=None, engine="python")
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    if selected_column not in df.columns:
        print(f"Error: column '{selected_column}' not found.")
        sys.exit(1)

    # Create categorized column
    df[result_column] = df[selected_column].apply(
        lambda x: "1 or below" if float(x) < selected_number else "2 or higher"
    )

    # Save output *as TSV*
    try:
        df.to_csv(output_file, sep="\t", index=False)
        print(f"TSV output written to {output_file}")
    except Exception as e:
        print(f"Error writing output: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()