#!/usr/bin/env python3

"""
Modified python script \
(09-02_extract_tss-feature_from_all-joined.assignedClusters_then_uniq.tss-feature.py)

Usage:
    python 09-02_uniq_tss_feature_modified.py \
    -i all-joined.assignedClusters.tsv
"""
import argparse
import pandas as pd

def main():
    """
    Extract TSS features from all-joined.assignedClusters.tsv
    """
    # コマンドライン引数の解析
    parser = argparse.ArgumentParser(description='TSS特徴抽出スクリプト')
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help='Path to input file (all-joined.assignedClusters.tsv)'
    )
    args = parser.parse_args()

    # 入力ファイルを読み込む
    input_file = args.input
    df = pd.read_csv(input_file, sep="\t", header=0)

    # 一時出力ファイル用のデータフレームを作成
    output_rows = []

    # 各行を処理
    for _, row in df.iterrows():
        # 1st, 2nd groupの 'tags.dominant_tss' を比較
        group_values = row[[7, 19]]  # 1st group: 7, 2nd group: 19

        # NAがあれば0に置き換え
        group_values = group_values.fillna(0)

        # 最大値を持つインデックスを取得
        max_value_idx = group_values.idxmax()

        # 最大値が0の場合はスキップ
        if group_values[max_value_idx] == 0:
            continue

        # インデックスを計算
        max_value_idxint = 0 if max_value_idx == 7 else 1

        # 最も高い 'tags.dominant_tss' を持つグループの情報を抽出
        group_chr = row.iloc[1 + max_value_idxint * 12]
        dominant_tss_value = row.iloc[5 + max_value_idxint * 12]
        group_cluster_start = row.iloc[2 + max_value_idxint * 12]
        group_cluster_end = row.iloc[3 + max_value_idxint * 12]
        group_tags = row.iloc[6 + max_value_idxint * 12]

        # dominant_tss が NaN の場合はスキップ
        if pd.isna(dominant_tss_value):
            continue

        dominant_tss = int(dominant_tss_value)  # dominant_tss を整数に変換
        group_cluster_start = int(group_cluster_start)  # cluster_start を整数に変換
        group_cluster_end = int(group_cluster_end)  # cluster_end を整数に変換

        # 'start'はdominant_tssから1を引いた位置に設定
        group_start = dominant_tss - 1  # 'start'はTSSの1塩基前

        # 最も高い 'tags.dominant_tss' を持つグループの 'strand', 'cluster', 'gene' を抽出
        group_strand = row.iloc[4 + max_value_idxint * 12]
        group_gene = row.iloc[11 + max_value_idxint * 12]

        # 出力行を作成
        output_row = [
            group_chr,
            'Cage_analysis', 
            'cage_cluster', 
            group_cluster_start,
            group_cluster_end, group_tags,
            'TSS', 
            group_start,
            dominant_tss,
            row.iloc[7 + max_value_idxint * 12],
            group_strand,
            row['cluster'],
            group_gene
        ]
        output_rows.append(output_row)

    # 一時データフレームを作成
    columns = [
        "chr",
        "source",
        "feature1",
        "start1", 
        "end1",
        "score1",
        "feature2",
        "start2",
        "end2",
        "tags.dominant_tss",
        "strand",
        "cluster",
        "gene"
    ]
    temp_df = pd.DataFrame(output_rows, columns=columns)

    # 'gene' 列を 'gene_id' に変更
    temp_df.rename(columns={'gene': 'gene_id'}, inplace=True)

    # 一時データフレームを "all_tss_feature.tsv" として保存
    temp_df.to_csv(
        "all_tss_feature.tsv",
        sep="\t",
        index=False
    )

    # 'tags.dominant_tss' を数値に変換（必要に応じて）
    temp_df['tags.dominant_tss'] = pd.to_numeric(temp_df['tags.dominant_tss'], errors='coerce')

    # 重複した遺伝子（gene_id）に対して、tags.dominant_tssが最も高い行を選択
    result_df = temp_df.loc[temp_df.groupby('gene_id')['tags.dominant_tss'].idxmax()]

    # 必要なカラムを選択して2つのデータフレームに分割
    df1 = result_df[['chr', 'source', 'feature1', 'start1', 'end1', 'score1', 'strand', 'cluster', 'gene_id']]
    df2 = result_df[['chr', 'source', 'feature2', 'start2', 'end2', 'tags.dominant_tss', 'strand', 'cluster', 'gene_id']]

    # 列のヘッダーを変更
    df1.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'tss_id', 'gene_id']
    df2.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'tss_id', 'gene_id']
    # 列の順番を変更し、end 列を2度使用
    df1 = df1[['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'end', 'source', 'score', 'gene_id']]
    # 6列目のヘッダーを「tss_pos」に変更
    df1.columns = ['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'tss_pos', 'source', 'score', 'gene_id']
    # 列の順番を変更し、end 列を2度使用
    df2 = df2[['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'end', 'source', 'score', 'gene_id']]
    # 6列目のヘッダーを「tss_pos」に変更
    df2.columns = ['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'tss_pos', 'source', 'score', 'gene_id']

    # 2度目に出現する 'end' 列を df2 の 'end' 列で置き換える
    df1['end_2'] = df2['end']  # 新しい列 'end_2' を作成して df2 の 'end' を代入

    # 列の順番を変更し、2度目の 'end' を df2 の値に置き換えたものにする
    df1 = df1[['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'end_2', 'source', 'score', 'gene_id']]

    # 6列目のヘッダーを「tss_pos」に変更
    df1.columns = ['tss_id', 'seqname', 'feature', 'start', 'end', 'strand', 'tss_pos', 'source', 'score', 'gene_id']

    # 1列目の値でソート
    df1 = df1.sort_values(by='tss_id')
    df2 = df2.sort_values(by='tss_id')

    # 分割したデータフレームをそれぞれ保存
    df1.to_csv("all_cage_cluster_feature_uniq.gene.tsv", sep="\t", index=False)
    df2.to_csv("all_tss_feature_uniq.gene.tsv", sep="\t", index=False)

    ### print("2つのデータフレームに分割して保存しました。")
    ### print("tss_positions.tsv と tss_metadata.tsv が生成されました。")

    ### # 列のヘッダーを変更
    ### result_df.columns = ['seqname', 'source', 'feature1', 'start1', 'end1', 'feature2', 'start2', 'end2', 'score', 'strand', 'tss_id', 'gene_id']
    ###
    ### # 列の順番を変更し、end 列を2度使用
    ### result_df = result_df[['tss_id', 'seqname', 'feature1', 'start1', 'end1', 'strand', 'end', 'source', 'feature', 'score', 'gene_id']]
    ###
    ### # 6列目のヘッダーを「tss_pos」に変更
    ### result_df.columns = ['tss_id', 'seqname', 'start', 'end', 'strand', 'tss_pos', 'source', 'feature', 'score', 'gene_id']
    ###
    ### # 1列目の値でソート
    ### result_df = result_df.sort_values(by='tss_id')
    ###
    ### # 結果を "all_tss_feature_uniq.gene.tsv" として保存
    ### result_df.to_csv("all_tss_feature_uniq.gene.tsv", sep="\t", index=False)
    ###
    ### print("処理が完了しました。")
    ### print("all_tss_feature.tsv と all_tss_feature_uniq.gene.tsv が生成されました。")

if __name__ == "__main__":
    main()
