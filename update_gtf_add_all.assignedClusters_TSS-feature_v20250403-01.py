import pandas as pd
import csv
import argparse
import os  # ファイル拡張子を取得するためのモジュール

def extract_gene_id(attribute_column, file_type):
    """
    GTFまたはGFFのattribute列からgene_idを抽出する関数
    """
    if file_type == "gtf":
        return attribute_column.str.extract(r'gene_id "([^"]+)"')
    elif file_type == "gff":
        # ID=gene-, Parent=gene-, または gene= の形式を探す
        return attribute_column.str.extract(r'(?:ID=gene-|Parent=gene-|gene=)([^;]+)')
    else:
        raise ValueError("Unsupported file type. Only GTF and GFF are supported.")

def extract_transcript_id(attribute_column, file_type):
    """
    GTFまたはGFFのattribute列からtranscript_idを抽出する関数
    """
    if file_type == "gtf":
        return attribute_column.str.extract(r'transcript_id "([^"]+)"')
    elif file_type == "gff":
        # Parent= の形式から transcript_id を抽出
        return attribute_column.str.extract(r'(?:ID=rna-|Parent=rna-|transcript_id=)([^;]+)')
    else:
        raise ValueError("Unsupported file type. Only GTF and GFF are supported.")

def update_gtf_with_tss(gtf_file, tss_file, output_file):
    # ファイルタイプを判定（拡張子で判断）
    file_extension = os.path.splitext(gtf_file)[1].lower()
    if file_extension == ".gtf":
        file_type = "gtf"
    elif file_extension == ".gff" or file_extension == ".gff3":
        file_type = "gff"
    else:
        raise ValueError("Unsupported file format. Please provide a .gtf or .gff file.")

    # ハッシュタグ行を保持するためにファイル全体を読み込む
    with open(gtf_file, 'r') as f:
        lines = f.readlines()

    # コメント行を抽出
    comment_lines = [(i, line) for i, line in enumerate(lines) if line.startswith('#')]

    # GTF/GFFファイルとTSSファイルの読み込み（コメント行を除外して読み込む）
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    tss_df = pd.read_csv(tss_file, sep='\t', header=0)

    # cDNA_match行を除外
    gtf_df = gtf_df[gtf_df['feature'] != 'cDNA_match']

    # gene_idとtranscript_idを抽出
    gtf_df['gene_id'] = extract_gene_id(gtf_df['attribute'], file_type)
    gtf_df['transcript_id'] = extract_transcript_id(gtf_df['attribute'], file_type)

    # デバッグ: gene_id の確認
    print("TSS gene IDs:", tss_df['gene_id'].unique())
    print("Total TSS gene IDs:", len(tss_df['gene_id'].unique()))
    print("GTF gene IDs:", gtf_df['gene_id'].unique())
    print("Total GTF gene IDs:", len(gtf_df['gene_id'].unique()))

    tss_gene_ids = tss_df['gene_id'].unique()
    updated_gtf = []  # 順番を保持するために list を使用

    # 全行のインデックスを保持
    gtf_df['original_index'] = gtf_df.index

    # gene_group に含まれなかった行を抽出
    processed_gene_ids = tss_gene_ids  # 処理対象の gene_id
    unprocessed_rows = gtf_df[~gtf_df['gene_id'].isin(processed_gene_ids)]

    # TSSのgene_idに基づいて更新
    for gene_id in tss_gene_ids:
        print(f"Processing gene_id: {gene_id}")
        gene_group = gtf_df[gtf_df['gene_id'] == gene_id]
        if gene_group.empty:
            print(f"Gene ID {gene_id} not found in GTF/GFF file.")
            continue

        cds = gene_group[gene_group['feature'] == 'CDS']
        if cds.empty:
            print(f"CDS is empty for gene_id: {gene_id}, but gene_group exists.")
            continue

        gene = gene_group[gene_group['feature'] == 'gene']
        exons = gene_group[gene_group['feature'] == 'exon']
        transcripts = gene_group[gene_group['feature'].isin(['transcript', 'mRNA'])]
        other_features = gene_group[~gene_group['feature'].isin(['gene', 'transcript', 'mRNA', 'exon', 'CDS'])]

        tss_group = tss_df[tss_df['gene_id'] == gene_id]

        if gene.empty or cds.empty:
            print(f"Skipping gene_id {gene_id} due to missing gene or CDS.")
            continue

        # gene_idに関連するtranscript_idを取得
        transcript_ids = cds['transcript_id'].dropna().unique()

        if len(transcript_ids) == 1:
            # 1つのtranscript_idを持つ場合の処理
            transcript_id = transcript_ids[0]
            print(f"Gene ID {gene_id} has a single transcript_id: {transcript_id}")
        else:
            # 複数のtranscript_idを持つ場合、TSSに最も近いtranscript_idを選択
            print(f"Gene ID {gene_id} has multiple transcript_ids: {transcript_ids}")
            tss_pos = tss_group['tss_pos'].iloc[0]

            # 各transcriptの開始点または終了点を取得し、TSSとの距離を計算
            transcript_distances = {}
            for transcript_id in transcript_ids:
                transcript_cds = cds[cds['transcript_id'] == transcript_id]
                if transcript_cds.empty:
                    continue
                if gene.iloc[0]['strand'] == '+':
                    transcript_start = transcript_cds['start'].min()
                else:
                    transcript_start = transcript_cds['end'].max()
                transcript_distances[transcript_id] = abs(transcript_start - tss_pos)

            # 最も近いtranscript_idを選択
            transcript_id = min(transcript_distances, key=transcript_distances.get)
            print(f"Selected transcript_id for gene_id {gene_id}: {transcript_id}")

        # 以下、選択されたtranscript_idを用いて処理を実行
        strand = gene.iloc[0]['strand']

        # gene の開始または終了位置を更新
        if strand == '+':
            new_gene_start = tss_group['tss_pos'].iloc[0]
            gene.iloc[0, gene.columns.get_loc('start')] = int(new_gene_start)
        else:
            new_gene_end = tss_group['tss_pos'].iloc[0]
            gene.iloc[0, gene.columns.get_loc('end')] = int(new_gene_end)

        updated_gtf.append({**gene.iloc[0].to_dict(), 'original_index': gene.index[0]})

        # exon の処理
        if strand == '+':
            first_exon_start = new_gene_start
        else:
            first_exon_end = new_gene_end

        if not exons.empty:
            first_exon = exons.sort_values(by='start' if strand == '+' else 'end', ascending=(strand == '+')).iloc[0].copy()
            if strand == '+':
                first_exon['start'] = first_exon_start
            else:
                first_exon['end'] = first_exon_end
            updated_gtf.append({**first_exon.to_dict(), 'original_index': first_exon.name})

            # 他の exon をそのまま追加
            other_exons = exons[~((exons['start'] == first_exon['start']) & (exons['end'] == first_exon['end']))]
            for idx, exon_row in other_exons.iterrows():
                updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})

        # transcript_id ごとに処理
        transcript_cds = cds[cds['transcript_id'] == transcript_id]

        # transcript_cds が空の場合はスキップ
        if transcript_cds.empty:
            print(f"No CDS found for transcript_id: {transcript_id}. Skipping.")
            continue

        if strand == '+':
            # 最初の CDS を基準に five_prime_utr_end を計算
            five_prime_utr_start = new_gene_start
            five_prime_utr_end = transcript_cds['start'].min() - 1
        else:
            # 最後の CDS を基準に five_prime_utr_start を計算
            five_prime_utr_start = transcript_cds['end'].max() + 1
            five_prime_utr_end = new_gene_end

        # five_prime_utr_start と five_prime_utr_end の範囲が正しい場合のみ追加
        if five_prime_utr_start <= five_prime_utr_end:
            # GTF形式かGFF形式かで処理を分岐
            if file_type == "gtf":
                # GTF形式: ダブルクォーテーションを付ける
                five_prime_utr_attribute = (
                    f'gene_id "{gene.iloc[0]["gene_id"]}"; '  # ダブルクォーテーション付き
                    f'Parent "rna-{transcript_id}"; '  # ダブルクォーテーション付き
                    f'five_prime_utr_id "tss_id-{tss_group.iloc[0]["tss_id"]}"; '  # ID=を追加
                    f'transcript_id "{transcript_id}"'  # ダブルクォーテーション付き
                )
            elif file_type == "gff":
                # GFF形式: ダブルクォーテーションを付けない
                five_prime_utr_attribute = (
                    f'ID=five_prime_utr-tss_id-{tss_group.iloc[0]["tss_id"]};'  # ID=を追加
                    f'gene_id={gene.iloc[0]["gene_id"]};'  # ダブルクォーテーションなし
                    f'Parent=rna-{transcript_id};'  # ダブルクォーテーションなし
                    f'transcript_id={transcript_id}'  # ダブルクォーテーションなし
                )
            else:
                raise ValueError("Unsupported file type. Only GTF and GFF are supported.")

            # 5'UTR行を作成
            five_prime_utr = pd.DataFrame({
                'seqname': [gene.iloc[0]['seqname']],
                'source': ['.'],
                'feature': ['five_prime_UTR'],
                'start': [int(five_prime_utr_start)],
                'end': [int(five_prime_utr_end)],
                'score': ['.'],
                'strand': [strand],
                'frame': ['.'],
                'attribute': [five_prime_utr_attribute]
            })

            # 5'UTR行をリストに追加
            updated_gtf.append({**five_prime_utr.iloc[0].to_dict(), 'original_index': gene.index[0]})
        else:
            print(f"Invalid UTR range for transcript_id: {transcript_id}. Skipping.")

        # transcript の処理（存在する場合のみ）
        if not transcripts.empty:
            updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
            if not updated_transcript.empty:
                if strand == '+':
                    updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(new_gene_start)
                else:
                    updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(new_gene_end)
                updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})

        # CDS の行をそのまま追加
        for idx, cds_row in cds.iterrows():
            updated_gtf.append({**cds_row.to_dict(), 'original_index': idx})

        # 他の要素をそのまま追加
        for idx, other_row in other_features.iterrows():
            # 10列目以降を削除
            other_row = other_row.iloc[:10]
            updated_gtf.append({**other_row.to_dict(), 'original_index': idx})

    # gene_group に含まれなかった行をそのまま追加
    for idx, row in unprocessed_rows.iterrows():
        updated_gtf.append({**row.to_dict(), 'original_index': idx})

    # 最終的に DataFrame に変換
    updated_gtf_df = pd.DataFrame(updated_gtf)

    # 元の順番を保持するためにソート
    updated_gtf_df = updated_gtf_df.sort_values(by='original_index').drop(columns=['original_index'])

    # 10列目以降を削除
    updated_gtf_df = updated_gtf_df.iloc[:, :9]

    # 重複を取り除く前にcDNA_match行を除外
    updated_gtf_df = updated_gtf_df[updated_gtf_df['feature'] != 'cDNA_match']
    updated_gtf_df = updated_gtf_df.drop_duplicates()

    # アップデートされた行と元の行を統合する処理を追加
    # 元のGTFデータと比較して、座標が変更された行を抽出
    # gtf_dfからcDNA_match行を除外してmerge処理を実行
    updated_gtf_only = updated_gtf_df.merge(
        gtf_df[gtf_df['feature'] != 'cDNA_match'].iloc[:, :9],  # cDNA_match行を除外
        how='left',
        indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    # 元の行の中で変更されなかった行を抽出
    unchanged_gtf = gtf_df[~gtf_df.index.isin(updated_gtf_only.index)]

    print("Unchanged GTF (CDS rows):")
    print(unchanged_gtf[unchanged_gtf['feature'] == 'CDS'])

    # アップデートされた行と変更されなかった行を結合
    final_gtf_df = pd.concat([updated_gtf_only, unchanged_gtf], ignore_index=True)

    # コメント行を元の位置に挿入
    output_lines = lines.copy()  # 元の行をコピー

    # 保存直前に元の順番を保持するためにソート
    final_gtf_df = final_gtf_df.sort_values(by='original_index').drop(columns=['original_index'])

    # コメント行を保持
    comment_lines = [line for line in output_lines if line.startswith('#')]

    # feature列のカスタム順序を定義
    feature_order = ['region', 'gene','transcript', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon', 'match', 'cDNA_match']
    final_gtf_df['feature'] = pd.Categorical(final_gtf_df['feature'], categories=feature_order, ordered=True)

    # 複数の列で一度にソート（featureをstartより優先）
    final_gtf_df = final_gtf_df.sort_values(
        by=['seqname', 'feature', 'start', 'end'],  # ソート条件
        ascending=[True, True, True, False]  # 昇順/降順の指定
    )

    final_gtf_df['gene_id'] = final_gtf_df['gene_id'].fillna('unknown_gene')
    final_gtf_df['gene_id'] = final_gtf_df['gene_id'].astype(str)

    # コメント行を先頭に追加し、データ部分を結合
    gtf_data_lines = final_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()
    output_lines = comment_lines + gtf_data_lines

    # 10列目以降を削除
    final_gtf_df = final_gtf_df.iloc[:, :9]  # DataFrameの9列目までを保持

    # データ部分を文字列形式に変換
    gtf_data_lines = final_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()

    # コメント行を先頭に追加し、データ部分を結合
    output_lines = comment_lines + gtf_data_lines

    # 空行を削除
    output_lines = [line.strip() for line in output_lines if line.strip()]

    # merge結果を確認
    merged_df = updated_gtf_df.merge(
        gtf_df.iloc[:, :9],  # 元のGTFデータ（座標情報のみ）
        how='left',
        indicator=True
    )
    print("Merge結果の内容:")
    print(merged_df['_merge'].value_counts())

    # left_onlyの行を確認
    print("Left onlyの行 (CDS rows):")
    print(merged_df[(merged_df['_merge'] == 'left_only') & (merged_df['feature'] == 'CDS')])

    # bothの行を確認
    print("Bothの行 (CDS rows):")
    print(merged_df[(merged_df['_merge'] == 'both') & (merged_df['feature'] == 'CDS')])

    print("Final GTF DataFrame (CDS rows):")
    print(final_gtf_df[final_gtf_df['feature'] == 'CDS'])

    # 保存処理
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')

###    cds_rows = final_gtf_df[final_gtf_df['feature'] == 'CDS']
###    cds_lines = cds_rows.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()
###    with open(output_file, 'a') as f:
###        f.write('\n'.join(cds_lines) + '\n')

if __name__ == "__main__":
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description="GTF/GFFファイルをTSS情報で更新するスクリプト")
    parser.add_argument("gtf_file", help="入力GTFまたはGFFファイルのパス")  # 修正
    parser.add_argument("tss_file", help="入力TSSファイルのパス")  # 修正
    parser.add_argument("output_file", help="出力GTFファイルのパス")  # 修正
    args = parser.parse_args()

    # 関数を実行
    update_gtf_with_tss(args.gtf_file, args.tss_file, args.output_file)