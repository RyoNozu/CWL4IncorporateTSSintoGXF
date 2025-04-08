import sys
import pandas as pd

def main():
    # Check if at least one file is provided
    if len(sys.argv) < 2:
        print("Usage: python join_all_assigned_clusters.py file1 file2 [file3 ... fileN]")
        print("Example: python join_all_assigned_clusters.py *.assignedClusters.txt")
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
        except Exception as e:
            print(f"Error processing file '{file}': {e}")
            sys.exit(1)

    # Merge all sorted dataframes
    try:
        merged_df = sorted_dataframes[0]
        for df in sorted_dataframes[1:]:
            merged_df = pd.merge(merged_df, df, on=merged_df.columns[0], how='outer', suffixes=('', '_dup'))
        merged_df.fillna("NA", inplace=True)
        print("All files merged successfully. First 5 rows of merged data:")
        print(merged_df.head())
    except Exception as e:
        print(f"Error during merging: {e}")
        sys.exit(1)

    # Save the merged dataframe to a file
    output_file = "all-joined.assignedClusters.tsv"
    try:
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"Merged data saved to '{output_file}'")
    except Exception as e:
        print(f"Error saving merged data to file: {e}")
        sys.exit(1)

# 使用例:
# このスクリプトは、複数のタブ区切りテキストファイル（例: *.assignedClusters.txt）を結合し、共通のキー列でマージします。
# 各ファイルは最初にキー列でソートされ、すべてのファイルが外部結合されます。
# 結果は "all-joined.assignedClusters.tsv" という名前のファイルに保存されます。
#
# 実行例:
# $ python 09-01_join_all_assignedClusters.py file1.assignedClusters.txt file2.assignedClusters.txt
# または
# $ python 09-01_join_all_assignedClusters.py *.assignedClusters.txt

if __name__ == "__main__":
    main()
