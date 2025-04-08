import pandas as pd

# 入力ファイルを読み込む
input_file = "all-joined.assignedClusters.tsv"
df = pd.read_csv(input_file, sep="\t", header=0)

# 一時出力ファイル用のデータフレームを作成
output_rows = []

# 各行を処理
for index, row in df.iterrows():
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

    # dominant_tss が NaN の場合はスキップ
    if pd.isna(dominant_tss_value):
        continue

    dominant_tss = int(dominant_tss_value)  # dominant_tss を整数に変換

    # 'start'はdominant_tssから1を引いた位置に設定
    group_start = dominant_tss - 1  # 'start'はTSSの1塩基前

    # 最も高い 'tags.dominant_tss' を持つグループの 'strand', 'cluster', 'gene' を抽出
    group_strand = row.iloc[4 + max_value_idxint * 12]
    group_gene = row.iloc[11 + max_value_idxint * 12]

    # 出力行を作成
    output_row = [
        group_chr, '.', 'TSS', group_start, dominant_tss, '.', group_strand, row['cluster'],
        row.iloc[7 + max_value_idxint * 12], group_gene
    ]
    output_rows.append(output_row)

# 一時データフレームを作成
columns = ["chr", "source", "feature", "start", "end", "score", "strand", "cluster", "tags.dominant_tss", "gene"]
temp_df = pd.DataFrame(output_rows, columns=columns)

# 'gene' 列を 'gene_id' に変更
temp_df.rename(columns={'gene': 'gene_id'}, inplace=True)

# 一時データフレームを "all_tss_feature.tsv" として保存
temp_df.to_csv("all_tss_feature.tsv", sep="\t", index=False)

# 'tags.dominant_tss' を数値に変換（必要に応じて）
temp_df['tags.dominant_tss'] = pd.to_numeric(temp_df['tags.dominant_tss'], errors='coerce')

# 重複した遺伝子（gene_id）に対して、tags.dominant_tssが最も高い行を選択
result_df = temp_df.loc[temp_df.groupby('gene_id')['tags.dominant_tss'].idxmax()]

# 列のヘッダーを変更
result_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'tss_id', 'tags.dominant_tss', 'gene_id']

# 列の順番を変更し、end 列を2度使用
result_df = result_df[['tss_id', 'seqname', 'start', 'end', 'strand', 'end', 'source', 'feature', 'score', 'tags.dominant_tss', 'gene_id']]

# 6列目のヘッダーを「tss_pos」に変更
result_df.columns = ['tss_id', 'seqname', 'start', 'end', 'strand', 'tss_pos', 'source', 'feature', 'score', 'tags.dominant_tss', 'gene_id']

# 1列目の値でソート
result_df = result_df.sort_values(by='tss_id')

# 結果を "all_tss_feature_uniq.gene.tsv" として保存
result_df.to_csv("all_tss_feature_uniq.gene.tsv", sep="\t", index=False)

print("処理が完了しました。")
print("all_tss_feature.tsv と all_tss_feature_uniq.gene.tsv が生成されました。")