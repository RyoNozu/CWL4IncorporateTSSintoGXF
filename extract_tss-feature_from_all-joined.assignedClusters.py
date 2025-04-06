#import pandas as pd
#
## 入力ファイルを読み込む
#input_file = "sort.all-joined.assignedClusters.tsv"
#df = pd.read_csv(input_file, sep="\t", header=0)
#
## 出力ファイルを開く
#output_file = "all_tss_feature.tsv"
#with open(output_file, "w") as out_file:
#    # ヘッダーを書き込む（もし必要であれば）
#    out_file.write("chr\tsource\tfeature\tstart\tend\tscore\tstrand\tcluster\ttags.dominant_tss\tgene\n")
#    
#    # 各行を処理
#    for index, row in df.iterrows():
#        # 各グループの 'tags.dominant_tss' を比較
#        group_values = row[7::12]  # 1st, 2nd, 3rd, 4th groupの 'tags.dominant_tss'
#
#        # NAがあれば0に置き換え
#        group_values = group_values.fillna(0)
#
#        # 最大値を持つインデックスを取得
#        max_value_idx = group_values.idxmax()
#
#        # 最大値が0の場合はスキップ
#        if group_values[max_value_idx] == 0:
#            continue
#
#        # 'tags.dominant_tss.2' から数字部分を抽出
#        max_value_idxint = int(max_value_idx.split('.')[-1])  # '2'を取り出して整数に変換
#
#        # 最も高い 'tags.dominant_tss' を持つグループの情報を抽出
#        group_chr = row[1 + max_value_idxint * 12]
#        dominant_tss = int(row[5 + max_value_idxint * 12])  # dominant_tss を整数に変換
#
#        # 'start'はdominant_tssから1を引いた位置に設定
#        group_start = dominant_tss - 1  # 'start'はTSSの1塩基前
#        group_end = int(row[3 + max_value_idxint * 12])  # 'end'列を整数に変換
#        group_strand = row[4 + max_value_idxint * 12]  # 'strand'列
#        group_max_value = row[7 + max_value_idxint * 12]  # 'tags.dominant_tss'列
#        group_gene = row[11 + max_value_idxint * 12]  # 'gene'列
#
#        # 出力内容を整形
#        output_line = f"{group_chr}\tTSSr\tTSS\t{group_start}\t{dominant_tss}\t.\t{group_strand}\t{row[0]}\t{group_max_value}\t{group_gene}\n"
#
#        # 結果をファイルに書き込む
#        out_file.write(output_line)


#import pandas as pd
#
## 入力ファイルを読み込む
#input_file = "sort.all-joined.assignedClusters.tsv"
#df = pd.read_csv(input_file, sep="\t", header=0)
#
## 出力ファイルを開く
#output_file = "all_tss_feature.tsv"
#with open(output_file, "w") as out_file:
#    # ヘッダーを書き込む（もし必要であれば）
#    out_file.write("chr\tsource\tfeature\tstart\tend\tscore\tstrand\tcluster\ttags.dominant_tss\tgene\n")
#    
#    # 各行を処理
#    for index, row in df.iterrows():
#        # 1st, 2nd groupの 'tags.dominant_tss' を比較
#        group_values = row[[7, 19]]  # 1st group: 7, 2nd group: 19
#
#        # NAがあれば0に置き換え
#        group_values = group_values.fillna(0)
#
#        # 最大値を持つインデックスを取得
#        max_value_idx = group_values.idxmax()
#
#        # 最大値が0の場合はスキップ
#        if group_values[max_value_idx] == 0:
#            continue
#
#        # インデックスを計算
#        max_value_idxint = 0 if max_value_idx == 7 else 1
#
#        # 最も高い 'tags.dominant_tss' を持つグループの情報を抽出
#        group_chr = row.iloc[1 + max_value_idxint * 12]
#        dominant_tss_value = row.iloc[5 + max_value_idxint * 12]
#
#        # dominant_tss が NaN の場合はスキップ
#        if pd.isna(dominant_tss_value):
#            continue
#
#        dominant_tss = int(dominant_tss_value)  # dominant_tss を整数に変換
#
#        # 'start'はdominant_tssから1を引いた位置に設定
#        group_start = dominant_tss - 1  # 'start'はTSSの1塩基前
#
#        # 出力行を作成
#        output_row = [
#            group_chr, '.', 'TSS', group_start, dominant_tss, '.', row['strand'], row['cluster'],
#            row.iloc[7 + max_value_idxint * 12], row['gene']
#        ]
#
#        # 出力ファイルに書き込む
#        out_file.write("\t".join(map(str, output_row)) + "\n")



import pandas as pd

# 入力ファイルを読み込む
input_file = "all-joined.assignedClusters.tsv"
df = pd.read_csv(input_file, sep="\t", header=0)

# 出力ファイルを開く
output_file = "all_tss_feature.tsv"
with open(output_file, "w") as out_file:
    # ヘッダーを書き込む（もし必要であれば）
    out_file.write("chr\tsource\tfeature\tstart\tend\tscore\tstrand\tcluster\ttags.dominant_tss\tgene\n")
    
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
        #group_cluster = row.iloc[0 + max_value_idxint * 12]
        group_gene = row.iloc[11 + max_value_idxint * 12]

        # 出力行を作成
        output_row = [
            group_chr, '.', 'TSS', group_start, dominant_tss, '.', group_strand, row['cluster'],
            row.iloc[7 + max_value_idxint * 12], group_gene
        ]

        # 出力ファイルに書き込む
        out_file.write("\t".join(map(str, output_row)) + "\n")