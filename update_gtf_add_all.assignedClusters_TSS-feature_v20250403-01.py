import pandas as pd
import csv
import argparse
import os  # ファイル拡張子を取得するためのモジュール
import re  # 正規表現モジュールをインポート

def extract_gene_id(attribute_column, file_type):
    """
    GTFまたはGFFのattribute列からgene_idを抽出する関数
    """
    if file_type == "gtf":
        return attribute_column.str.extract(r'gene_id\s*"?([^";]+)"?')
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
    if (file_extension == ".gtf"):
        file_type = "gtf"
    elif (file_extension == ".gff") or (file_extension == ".gff3"):
        file_type = "gff"
    else:
        raise ValueError("Unsupported file format. Please provide a .gtf or .gff file.")

    # ハッシュタグ行を保持するためにファイル全体を読み込む
    with open(gtf_file, 'r') as f:
        lines = f.readlines()

    # コメント行を抽出（文字列部分のみ保持）
    comment_lines = [line.strip() for line in lines if line.startswith('#') and line.strip() != '']

    # GTF/GFFファイルとTSSファイルの読み込み（コメント行を除外して読み込む）
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    tss_df = pd.read_csv(tss_file, sep='\t', header=0)

    # cDNA_match行を除外
    gtf_df = gtf_df[gtf_df['feature'] != 'cDNA_match']

    # gene_idとtranscript_idを抽出
    gtf_df['gene_id'] = extract_gene_id(gtf_df['attribute'], file_type)
    print(gtf_df[['gene_id', 'attribute']].head(10))
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
        transcript_ids = transcripts['transcript_id'].dropna().unique()

        if len(transcript_ids) == 1:
            # 1つのtranscript_idを持つ場合の処理
            transcript_id = transcript_ids[0]
            print(f"Gene ID {gene_id} has a single transcript_id: {transcript_id}")

            # TSS位置を取得
            tss_pos = tss_group['tss_pos'].iloc[0]
            strand = gene.iloc[0]['strand']

            # gene start または gene end の更新
            if strand == '+':
                new_gene_start = int(tss_pos)
                gene.iloc[0, gene.columns.get_loc('start')] = new_gene_start
                print(f"Updated gene start for gene_id {gene_id} to {tss_pos}.")
            else:
                new_gene_end = int(tss_pos)
                gene.iloc[0, gene.columns.get_loc('end')] = new_gene_end
                print(f"Updated gene end for gene_id {gene_id} to {tss_pos}.")

            # transcript の更新処理
            if not transcripts.empty:
                updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                if not updated_transcript.empty:
                    if strand == '+':
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(tss_pos)
                    else:
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(tss_pos)
                    updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                    print(f"Updated transcript for transcript_id {transcript_id}:")
                    print(updated_transcript.iloc[0])

            # exon の更新処理
            if not exons.empty:
                related_exons = exons[exons['attribute'].str.contains(f'transcript_id={transcript_id}')]

                if strand == '+':
                    # start が最小の exon を取得
                    smallest_start_exon = related_exons.loc[related_exons['start'].idxmin()].copy()
                    smallest_start_exon['start'] = int(tss_pos)  # TSS に基づいて start を更新
                    updated_gtf.append({**smallest_start_exon.to_dict(), 'original_index': smallest_start_exon.name})
                    print(f"Updated exon for transcript_id {transcript_id} (smallest start):")
                    print(smallest_start_exon)

                    # 他の exon をそのまま追加
                    other_exons = related_exons[related_exons.index != smallest_start_exon.name]
                    for idx, exon_row in other_exons.iterrows():
                        updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                        print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        print(exon_row)

                else:
                    # end が最大の exon を取得
                    largest_end_exon = related_exons.loc[related_exons['end'].idxmax()].copy()
                    largest_end_exon['end'] = int(tss_pos)  # TSS に基づいて end を更新
                    updated_gtf.append({**largest_end_exon.to_dict(), 'original_index': largest_end_exon.name})
                    print(f"Updated exon for transcript_id {transcript_id} (largest end):")
                    print(largest_end_exon)

                    # 他の exon をそのまま追加
                    other_exons = related_exons[related_exons.index != largest_end_exon.name]
                    for idx, exon_row in other_exons.iterrows():
                        updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                        print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        print(exon_row)

        else:
            # 複数のtranscript_idを持つ場合、TSSに最も近いtranscript_idを選択
            print(f"Gene ID {gene_id} has multiple transcript_ids: {transcript_ids}")
            tss_pos = tss_group['tss_pos'].iloc[0]

            # 各transcriptのCDS開始点または終了点を取得し、TSSとの距離を計算
            tss_to_cds_distance = {}  # 修正: transcript_distances → tss_to_cds_distance
            for transcript_id in transcript_ids:
                transcript_cds = cds[cds['transcript_id'] == transcript_id]
                if transcript_cds.empty:
                    continue
                if gene.iloc[0]['strand'] == '+':
                    five_utr_end = transcript_cds['start'].min()  # 修正: transcript_start → five_utr_end
                else:
                    five_utr_end = transcript_cds['end'].max()  # 修正: transcript_start → five_utr_end
                tss_to_cds_distance[transcript_id] = abs(five_utr_end - tss_pos)  # 修正: transcript_distances → tss_to_cds_distance

            # 最も近い transcript_id を選択
            transcript_id = min(tss_to_cds_distance, key=tss_to_cds_distance.get)
            print(f"Selected transcript_id for gene_id {gene_id}: {transcript_id}")

        # 選ばれた transcript_id の情報を取得
        transcript_cds = cds[cds['transcript_id'] == transcript_id]
        transcript_start = transcripts['start'].min() if not transcripts.empty else None

        # 以下、選択されたtranscript_idを用いて処理を実行
        strand = gene.iloc[0]['strand']

        # gene の開始または終了位置を更新
        if strand == '+':
            # TSSに最も近いtranscriptのstartを取得
            tss_pos = tss_group['tss_pos'].iloc[0]

            # 各transcriptの開始位置とTSSの距離を計算
            transcript_distances = {}
            for idx, transcript_row in transcripts.iterrows():
                transcript_start = transcript_row['start']
                transcript_distances[transcript_start] = abs(transcript_start - tss_pos)

            # デバッグ: transcript_distances を確認
            print(f"Debug: transcript_distances for gene_id {gene_id}: {transcript_distances}")

            # TSSに最も近いtranscriptのstartを取得
            closest_transcript_start = min(transcript_distances, key=transcript_distances.get)

            # gene の開始位置を更新
            print(f"Debug: closest_transcript_start = {closest_transcript_start}, gene_start = {gene.iloc[0]['start']}")
            if closest_transcript_start == gene.iloc[0]['start']:
                new_gene_start = int(tss_pos)
                gene.iloc[0, gene.columns.get_loc('start')] = new_gene_start
                print(f"Updated gene start for gene_id {gene_id} to {tss_pos} (matched closest transcript start).")
            else:
                new_gene_start = gene.iloc[0]['start']
                print(f"Gene start for gene_id {gene_id} not updated (closest transcript start does not match gene start).")
            
            # transcript の更新処理を追加
            if not transcripts.empty:
                updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                if not updated_transcript.empty:
                    if strand == '+':
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(tss_pos)
                    else:
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(tss_pos)
                    updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                    print(f"Updated transcript for transcript_id {transcript_id}:")
                    print(updated_transcript.iloc[0])

                # 選択されなかった transcript_id の行をそのまま追加
                other_transcripts = transcripts[transcripts['transcript_id'] != transcript_id]
                for idx, transcript_row in other_transcripts.iterrows():
                    updated_gtf.append({**transcript_row.to_dict(), 'original_index': idx})
                    print(f"Unmodified transcript added for transcript_id {transcript_row['transcript_id']}:")
                    print(transcript_row)

            # exon の更新処理を追加
            if not exons.empty:
                # 選択された transcript_id に関連する exon を取得
                related_exons = exons[exons['attribute'].str.contains(f'transcript_id={transcript_id}')]

                if strand == '+':
                    # start が最小の exon を取得
                    smallest_start_exon = related_exons.loc[related_exons['start'].idxmin()].copy()
                    smallest_start_exon['start'] = int(tss_pos)  # TSS に基づいて start を更新
                    updated_gtf.append({**smallest_start_exon.to_dict(), 'original_index': smallest_start_exon.name})
                    print(f"Updated exon for transcript_id {transcript_id} (smallest start):")
                    print(smallest_start_exon)

                    # 他の exon をそのまま追加
                    other_exons = related_exons[related_exons.index != smallest_start_exon.name]
                    for idx, exon_row in other_exons.iterrows():
                        updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                        print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        print(exon_row)

                else:
                    # strand が '-' の場合、end が最大の exon を取得
                    largest_end_exon = related_exons.loc[related_exons['end'].idxmax()].copy()
                    largest_end_exon['end'] = int(tss_pos)  # TSS に基づいて end を更新

                    # 更新された exon を追加
                    updated_gtf.append({**largest_end_exon.to_dict(), 'original_index': largest_end_exon.name})
                    print(f"Updated exon for transcript_id {transcript_id} (largest end):")
                    print(largest_end_exon)

                    # 他の exon をそのまま追加
                    other_exons = related_exons[related_exons.index != largest_end_exon.name]
                    for idx, exon_row in other_exons.iterrows():
                        updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                        print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        print(exon_row)

        else:
            # TSSに最も近いtranscriptのendを取得
            transcript_end = transcript_cds['end'].max()  
            tss_pos = tss_group['tss_pos'].iloc[0]

            # TSSに最も近いtranscriptのendがgene_endと一致する場合のみgene_endを更新
            if  transcript_end == gene.iloc[0]['end']:  
                new_gene_end = int(tss_pos)  # 更新されたgene_endをnew_gene_endに設定
                gene.iloc[0, gene.columns.get_loc('end')] = new_gene_end
                print(f"Updated gene end for gene_id {gene_id} to {tss_pos} (matched transcript end).")
            else:
                new_gene_end = gene.iloc[0]['end']  # 元のgene_endを使用
                print(f"Gene end for gene_id {gene_id} not updated (transcript end does not match gene end).")

            # transcript の更新処理を追加
            if not transcripts.empty:
                updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                if not updated_transcript.empty:
                    if strand == '+':
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(new_gene_start)
                    else:
                        updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(new_gene_end)
                    updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                    print(f"Updated transcript for transcript_id {transcript_id}:")
                    print(updated_transcript.iloc[0])

                # 選択されなかった transcript_id の行をそのまま追加
                other_transcripts = transcripts[transcripts['transcript_id'] != transcript_id]
                for idx, transcript_row in other_transcripts.iterrows():
                    updated_gtf.append({**transcript_row.to_dict(), 'original_index': idx})
                    print(f"Unmodified transcript added for transcript_id {transcript_row['transcript_id']}:")
                    print(transcript_row)

            # exon の更新処理を追加
            if not exons.empty:
                # 選択された transcript_id に関連する exon を取得
                related_exons = exons[exons['attribute'].str.contains(f'transcript_id={transcript_id}')]

                for idx, exon_row in related_exons.iterrows():
                    updated_exon = exon_row.copy()
                    if strand == '+':
                        updated_exon['start'] = int(tss_pos)  # TSS に基づいて start を更新
                    else:
                        updated_exon['end'] = int(tss_pos)  # TSS に基づいて end を更新

                    # 更新された exon を追加
                    updated_gtf.append({**updated_exon.to_dict(), 'original_index': idx})
                    print(f"Updated exon for transcript_id {transcript_id}:")
                    print(updated_exon)

                # 選択されなかった exon をそのまま追加
                other_exons = exons[~exons['attribute'].str.contains(f'transcript_id={transcript_id}')]
                for idx, exon_row in other_exons.iterrows():
                    updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                    print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                    print(exon_row)

        updated_gtf.append({**gene.iloc[0].to_dict(), 'original_index': gene.index[0]})

        if strand == '+':
            # 最初の CDS を基準に five_prime_utr_end を計算
            five_prime_utr_start = tss_pos  # 修正: new_gene_start → tss_pos
            five_prime_utr_end = transcript_cds['start'].min() - 1
        else:
            # 最後の CDS を基準に five_prime_utr_start を計算
            five_prime_utr_start = transcript_cds['end'].max() + 1 
            five_prime_utr_end = tss_pos  # 修正: new_gene_end → tss_pos

        # five_prime_UTR_start と five_prime_UTR_end の範囲が正しい場合のみ追加
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

            # 対応する transcript_id の最後の CDS 行のインデックスを取得
            last_cds_index = cds[cds['transcript_id'] == transcript_id].index.max()

            # five_prime_UTR行をリストに追加
            five_prime_utr_dict = {**five_prime_utr.iloc[0].to_dict(), 'original_index': last_cds_index}
            updated_gtf.append(five_prime_utr_dict)

            # デバッグ: five_prime_UTR の情報を確認
            print("Added five_prime_UTR entry:")
            print(five_prime_utr_dict)
        else:
            print(f"Invalid UTR range for transcript_id: {transcript_id}. Skipping.")

        # transcript の処理（存在する場合のみ）
        if not transcripts.empty:
            updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
            if not updated_transcript.empty:  # 修正: () を削除
                if strand == '+':
                    updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(tss_pos)
                else:
                    updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(tss_pos)
                
                # デバッグ: 更新された transcript の内容を確認
                print(f"Debug: Updated transcript for transcript_id {transcript_id}:")
                print(updated_transcript.iloc[0])
                
                updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                
                # デバッグ: updated_gtf に追加された内容を確認
                print(f"Debug: Updated GTF entry for transcript_id {transcript_id}:")
                print(updated_gtf[-1])

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

    # デバッグ: updated_gtf に five_prime_UTR が含まれているか確認
    print("Debug: Updated GTF entries (first 10):")
    print(updated_gtf[:50])  # 最初の10行を確認

    # 最終的に DataFrame に変換
    updated_gtf_df = pd.DataFrame(updated_gtf)

    # デバッグ: updated_gtf_df に five_prime_UTR が含まれているか確認
    print("Debug: Updated GTF DataFrame preview:")
    print(updated_gtf_df[updated_gtf_df['feature'] == 'five_prime_UTR'])

    # 重複を取り除く前にcDNA_match行を除外
    updated_gtf_only = updated_gtf_df.drop_duplicates()

    # デバッグ: updated_gtf_only に five_prime_UTR が含まれているか確認
    print("Debug: Updated GTF only (five_prime_UTR):")
    print(updated_gtf_only[updated_gtf_only['feature'] == 'five_prime_UTR'])

    # 元の行の中で変更されなかった行を抽出
    unchanged_gtf = gtf_df[~gtf_df.index.isin(updated_gtf_only.index)]

    # アップデートされた行と変更されなかった行を結合
    final_gtf_df = pd.concat([updated_gtf_only, unchanged_gtf], ignore_index=True)

    # デバッグ: final_gtf_df に five_prime_UTR が含まれているか確認
    print("Debug: Final GTF DataFrame (five_prime_UTR):")
    print(final_gtf_df[final_gtf_df['feature'] == 'five_prime_UTR'])

    # コメント行を先頭に追加し、データ部分を結合する前にソート処理を実行
    # feature列のカスタム順序を設定
    feature_order = ['region', 'gene', 'transcript', 'mRNA', 'exon', 'CDS', 
                     'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon', 
                     'match', 'cDNA_match']
    final_gtf_df['feature'] = pd.Categorical(final_gtf_df['feature'], categories=feature_order, ordered=True)

    # original_index と feature でソート
    final_gtf_df = final_gtf_df.sort_values(by=['original_index', 'feature'], ascending=[True, True])

    # データ部分を文字列形式に変換
    gtf_data_lines = final_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()

    # コメント行を先頭に追加し、データ部分を結合
    output_lines = comment_lines + gtf_data_lines

    # デバッグ: output_lines に five_prime_UTR が含まれているか確認
    print("Debug: Output lines (five_prime_UTR):")
    print([line for line in output_lines if 'five_prime_UTR' in line])

    # 保存処理
    # 正規表現を使用して空行を削除
    output_lines = [line for line in output_lines if not re.match(r'^\s*$', line)]

    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')

    print("File saved successfully.")

    print(gtf_df[['attribute']].head(50))
    print("Debug: Final GTF DataFrame preview:")
    print(final_gtf_df.head(50))  # 最初の10行を確認

if __name__ == "__main__":
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description="GTF/GFFファイルをTSS情報で更新するスクリプト")
    parser.add_argument("gtf_file", help="入力GTFまたはGFFファイルのパス")  # 修正
    parser.add_argument("tss_file", help="入力TSSファイルのパス")  # 修正
    parser.add_argument("output_file", help="出力GTFファイルのパス")  # 修正
    args = parser.parse_args()

    # 関数を実行
    update_gtf_with_tss(args.gtf_file, args.tss_file, args.output_file)