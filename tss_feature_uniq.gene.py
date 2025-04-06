import pandas as pd

df = pd.read_csv("all_tss_feature.tsv", sep='\t')

# 'tags.dominant_tss' を数値に変換（必要に応じて）
df['tags.dominant_tss'] = pd.to_numeric(df['tags.dominant_tss'], errors='coerce')

# 重複した遺伝子（gene）に対して、tags.dominant_tssが最も高い行を選択
result_df = df.loc[df.groupby('gene')['tags.dominant_tss'].idxmax()]

# 列のヘッダーを変更
result_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'tss_id', 'tags.dominant_tss', 'gene']

# 列の順番を変更し、end 列を2度使用
result_df = result_df[['tss_id', 'seqname', 'start', 'end', 'strand', 'end', 'source', 'feature', 'score', 'tags.dominant_tss', 'gene']]

# 6列目のヘッダーを「tss_pos」に変更
result_df.rename(columns={'end': 'tss_pos'}, inplace=True)

# 1列目の値でソート
result_df = result_df.sort_values(by='tss_id')

# 結果をタブ区切りのテキストファイルとして保存
output_file = 'all_tss_feature_uniq.gene.tsv'
result_df.to_csv(output_file, sep='\t', index=False)

print(f"処理が完了しました。結果は {output_file} に保存されました。")
