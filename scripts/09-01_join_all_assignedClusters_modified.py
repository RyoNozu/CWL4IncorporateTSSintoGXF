#! /usr/bin/env python3

"""
Usage:
This script joins multiple tab-separated text files \
    (e.g., *.assignedClusters.txt) by a common key column.
Each file is sorted by the key column, and all files are merged using an outer join.
The result is saved to a file named "all-joined.assignedClusters.tsv".

Example:
$ python3 09-01_join_all_assignedClusters.py file1.assignedClusters.txt file2.assignedClusters.txt
or
$ python3 09-01_join_all_assignedClusters.py *.assignedClusters.txt
"""
import sys
import pandas as pd

def main():
    """
    This function joins multiple tab-separated text files \
        (e.g., *.assignedClusters.txt) by a common key column.
    Each file is sorted by the key column, and all files are merged using an outer join.
    The result is saved to a file named "all-joined.assignedClusters.tsv".
    """
    # Check if at least one file is provided
    if len(sys.argv) < 2:
        print("Usage: python3 09-01_join_all_assignedClusters.py file1 file2 [file3 ... fileN]")
        print("Example: python3 09-01_join_all_assignedClusters.py *.assignedClusters.txt")
        sys.exit(1)

    # List of input files
    input_files = sys.argv[1:]
    print(f"Input files: {input_files}")

    # Read and sort each file
    sorted_dataframes = []
    for file in input_files:
        try:
            # Read the file
            df = pd.read_csv(file, sep="\t", header=0)
            print(f"File '{file}' read successfully. First 5 rows:")
            print(df.head())

            # Ensure the first column (key column) is of consistent type (e.g., string)
            key_column = df.columns[0]
            df[key_column] = df[key_column].astype(str)

            # Sort by the first column (assumed to be the key column)
            df_sorted = df.sort_values(by=key_column)
            sorted_dataframes.append(df_sorted)
            print(f"File '{file}' sorted successfully. First 5 rows of sorted data:")
            print(df_sorted.head())
        except (pd.errors.ParserError, IOError, ValueError, KeyError) as e:
            print(f"Error processing file '{file}': {e}")
            sys.exit(1)

    # Merge all sorted dataframes
    try:
        merged_df = sorted_dataframes[0]
        for df in sorted_dataframes[1:]:
            merged_df = pd.merge(
                merged_df, df, on=merged_df.columns[0], how='outer', suffixes=('', '_dup')
            )
        merged_df.fillna("NA", inplace=True)
        print("All files merged successfully. First 5 rows of merged data:")
        print(merged_df.head())
    except (ValueError, KeyError, IndexError) as e: # if the files are not mergeable, raise an error
        print(f"Error during merging: {e}")
        sys.exit(1)

    # Save the merged dataframe to a file
    output_file = "all-joined.assignedClusters.tsv"
    try:
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"Merged data saved to '{output_file}'")
    except IOError as e: # if the file is not writable
        print(f"Error saving merged data to file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
