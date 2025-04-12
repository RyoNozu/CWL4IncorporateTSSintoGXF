import pandas as pd
import csv
import argparse
import os  # ファイル拡張子を取得するためのモジュール
import re  # 正規表現モジュールをインポート

# 出力の幅を広げる
###pd.set_option('display.max_colwidth', None)  # 列の最大幅を無制限に設定
pd.set_option('display.max_columns', 250)  # 表示する列数を無制限に設定
pd.set_option('display.width', 250)        # 出力の幅を広げる
pd.set_option('display.max_rows', 250)     # 表示する行数を無制限に設定

def extract_gene_id(attribute_column, file_type):
    """
    GTFまたはGFFのattribute列からgene_idを抽出する関数
    """
    if file_type == "gtf":
        return attribute_column.str.extract(r'gene_id\s*"?([^";]+)"?')
    elif file_type == "gff":
        # ID=gene-, Parent=gene-, または gene= の形式を探す
        return attribute_column.str.extract(r'(?:gene_id=|ID=gene-|Parent=gene-|gene=)([^;]+)')
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
        return attribute_column.str.extract(r'(?:transcript_id=|ID=rna-|Parent=rna-)([^;]+)')
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

    # オリジナルのgtf_dfのスタッツを確認するためにgtf_dfのコピー
    gtf_df_original = gtf_df.copy()  # 元のデータフレームを保持

    # gene_idとtranscript_idを抽出
    gtf_df['gene_id'] = extract_gene_id(gtf_df['attribute'], file_type)
    print(gtf_df[['gene_id', 'attribute']].head(75))
    gtf_df['transcript_id'] = extract_transcript_id(gtf_df['attribute'], file_type)
    print(gtf_df[['transcript_id', 'attribute']].head(75))

    # デバッグ: gene_id の確認
    print("TSS gene IDs:", tss_df['gene_id'].unique())
    print("Total TSS gene IDs:", len(tss_df['gene_id'].unique()))
    print("GTF gene IDs:", gtf_df['gene_id'].unique())
    print("Total GTF gene IDs:", len(gtf_df['gene_id'].unique()))

    tss_gene_ids = tss_df['gene_id'].unique()
    updated_gtf = []  # 順番を保持するために list を使用
    unchanged_gtf = []  # 更新されなかった行を保持するためのリスト
    cds_gtf = []  # CDS行を保持するためのリスト

    # 全行のインデックスを保持
    gtf_df['original_index'] = gtf_df.index

    # gene_group に含まれなかった行を抽出
    processed_gene_ids = tss_gene_ids  # 処理対象の gene_id
    unprocessed_rows = gtf_df[
        ~gtf_df['gene_id'].isin(processed_gene_ids)
    ]

    # デバッグ: unprocessed_rows の内容を出力
    print("\nDebug: Unprocessed rows (rows with gene_id not in processed_gene_ids):")
    print(unprocessed_rows)
    
    # 必要に応じてファイルに出力
    unprocessed_rows_output_file = "unprocessed_rows_debug_output.tsv"
    unprocessed_rows.to_csv(unprocessed_rows_output_file, sep='\t', index=False)
    print(f"Unprocessed rows have been written to {unprocessed_rows_output_file}")

    # 更新の対象とならない行を unchanged_gtf に追加
    for idx, row in unprocessed_rows.iterrows():
        unchanged_gtf.append({**row.to_dict(), 'original_index': idx})

    # TSSのgene_idに基づいて更新
    for gene_id in tss_gene_ids:
        print(f"Processing gene_id: {gene_id}")
        gene_group = gtf_df[gtf_df['gene_id'] == gene_id]
        if gene_group.empty:
            print(f"Gene ID {gene_id} not found in GTF/GFF file.")
            continue

        cds = gene_group[gene_group['feature'] == 'CDS']
        gene = gene_group[gene_group['feature'].isin(['gene', 'pseudogene'])]

        if gene.empty or cds.empty:
            print(f"Skipping gene_id {gene_id} due to missing gene or CDS.")
            
            # 該当する gene_id を持つ行を unchanged_gtf に追加
            for idx, row in gene_group.iterrows():
                unchanged_gtf.append({**row.to_dict(), 'original_index': idx})
                print(f"Added row for gene_id {gene_id} to unchanged_gtf: {row.to_dict()}")
            
            continue

        # CDS の行をそのままcds_gtfに追加
        for idx, cds_row in cds.iterrows():
            cds_gtf.append({**cds_row.to_dict(), 'original_index': idx})

        exons = gene_group[gene_group['feature'] == 'exon']
        transcripts = gene_group[gene_group['feature'].isin(['transcript', 'mRNA', 'lnc_RNA', 'miRNA', 'tRNA', 'rRNA', 'antisense_RNA'])]
        other_features = gene_group[~gene_group['feature'].isin(['gene', 'pseudogene', 'transcript', 'mRNA', 'exon', 'CDS', 'lnc_RNA', 'miRNA', 'tRNA', 'rRNA', 'antisense_RNA'])]

        tss_group = tss_df[tss_df['gene_id'] == gene_id]

        # 更新の対象となる可能性の無い行(feature)を unchanged_gtf に追加
        for idx, row in other_features.iterrows():
            unchanged_gtf.append({**row.to_dict(), 'original_index': idx})

        # 各gene_idと紐づいたtranscript_idをexon行から取得(transcript, CDS, exon 行いずれもOK)
        transcript_ids = exons['transcript_id'].dropna().unique()
        print(f"Transcript IDs for gene_id {gene_id}: {transcript_ids}")

        if len(transcript_ids) == 1:
            # 1つのtranscript_idを持つ場合の処理
            transcript_id = transcript_ids[0]
            print(f"Gene ID {gene_id} has a single transcript_id: {transcript_id}")

            # TSS位置を取得
            tss_pos = tss_group['tss_pos'].iloc[0]
            strand = gene.iloc[0]['strand']

            # strandに基づいた処理
            if strand == '+':
                if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
                    # gene 座標の更新
                    new_gene_start = int(tss_pos)
                    gene.iloc[0, gene.columns.get_loc('start')] = new_gene_start
                    print(f"Updated gene start for gene_id {gene_id} to {tss_pos}.")
                    # 更新された gene 行を updated_gtf に追加
                    updated_gene = gene.iloc[0].to_dict()
                    updated_gene['original_index'] = gene.index[0]
                    updated_gtf.append(updated_gene)
                else:
                    print(f"Gene start for gene_id {gene_id} not updated (TSS position is inside CDS).")
                    # オリジナルの gene 行を updated_gtf に追加
                    original_gene = gene.iloc[0].to_dict()
                    original_gene['original_index'] = gene.index[0]
                    unchanged_gtf.append(original_gene)
                    print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")

                # transcript の更新処理
                if not transcripts.empty:
                    updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                    if not updated_transcript.empty:
                        if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
                            updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(tss_pos)
                            updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                            print(f"Updated transcript for transcript_id {transcript_id}:")
                            print(updated_transcript.iloc[0])
                        else:
                            print(f"Transcript start for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                            # オリジナルの transcript or mRNA 行を updated_gtf に追加
                            updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].iloc[0].to_dict()
                            updated_transcript['original_index'] = transcripts[transcripts['transcript_id'] == transcript_id].index[0]
                            unchanged_gtf.append(updated_transcript)
                            print(f"Original transcript entry added to updated_gtf for transcript_id {transcript_id}:")

                # exon の更新処理
                related_exons = exons[exons['attribute'].str.contains(
                    rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                )]
                if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
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
                        ###print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        ###print(exon_row)
                else:
                    print(f"Exon start for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                    # オリジナルの exon 行を updated_gtf に追加
                    original_exon = related_exons.iloc[0].to_dict()
                    original_exon['original_index'] = related_exons.index[0]
                    unchanged_gtf.append(original_exon)
                    print(f"Original exon entry added to updated_gtf for transcript_id {transcript_id}:")

                # five_prime_UTR の計算と追加処理
                # 最初の CDS を基準に five_prime_utr_end を計算
                five_prime_utr_start = tss_pos
                five_prime_utr_end = cds[cds['transcript_id'] == transcript_id]['start'].min() - 1

                # five_prime_UTR_start と five_prime_UTR_end の範囲が正しい場合のみ追加
                if five_prime_utr_start <= five_prime_utr_end:
                    if file_type == "gtf":
                        five_prime_utr_attribute = (
                            f'gene_id "{gene.iloc[0]["gene_id"]}"; '
                            f'Parent "rna-{transcript_id}"; '
                            f'five_prime_utr_id "tss_id-{tss_group.iloc[0]["tss_id"]}"; '
                            f'transcript_id "{transcript_id}"'
                        )
                    elif file_type == "gff":
                        five_prime_utr_attribute = (
                            f'ID=five_prime_utr-tss_id-{tss_group.iloc[0]["tss_id"]};'
                            f'gene_id={gene.iloc[0]["gene_id"]};'
                            f'Parent=rna-{transcript_id};'
                            f'transcript_id={transcript_id}'
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

            # strand: '-' の処理
            else:
                if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
                    new_gene_end = int(tss_pos)
                    gene.iloc[0, gene.columns.get_loc('end')] = new_gene_end
                    print(f"Updated gene end for gene_id {gene_id} to {tss_pos}.")
                    # 更新された gene 行を updated_gtf に追加
                    updated_gene = gene.iloc[0].to_dict()
                    updated_gene['original_index'] = gene.index[0]
                    updated_gtf.append(updated_gene)
                else:
                    print(f"Gene start for gene_id {gene_id} not updated (TSS position is inside CDS).")
                    # オリジナルの gene 行を updated_gtf に追加
                    original_gene = gene.iloc[0].to_dict()
                    original_gene['original_index'] = gene.index[0]
                    unchanged_gtf.append(original_gene)
                    print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")
                
                # transcript の更新処理
                if not transcripts.empty:
                    updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                    if not updated_transcript.empty:
                        if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
                            updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(tss_pos)
                            updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                            print(f"Updated transcript for transcript_id {transcript_id}:")
                            print(updated_transcript.iloc[0])
                        else:
                            print(f"Transcript end for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                            # オリジナルの transcript or mRNA 行を updated_gtf に追加
                            updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].iloc[0].to_dict()
                            updated_transcript['original_index'] = transcripts[transcripts['transcript_id'] == transcript_id].index[0]
                            unchanged_gtf.append(updated_transcript)
                            print(f"Original transcript entry added to updated_gtf for transcript_id {transcript_id}:")

                # exon の更新処理
                related_exons = exons[exons['attribute'].str.contains(
                    rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                )]
                if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
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
                        ###print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                        ###print(exon_row)
                else:
                    print(f"Exon end for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                    # オリジナルの exon 行を updated_gtf に追加
                    original_exon = related_exons.iloc[0].to_dict()
                    original_exon['original_index'] = related_exons.index[0]
                    unchanged_gtf.append(original_exon)
                    print(f"Original exon entry added to updated_gtf for transcript_id {transcript_id}:")

                # five_prime_UTR の計算と追加処理
                # 最後の CDS を基準に five_prime_utr_start を計算
                five_prime_utr_start = cds[cds['transcript_id'] == transcript_id]['end'].max() + 1
                five_prime_utr_end = tss_pos

                # five_prime_UTR_start と five_prime_UTR_end の範囲が正しい場合のみ追加
                if five_prime_utr_start <= five_prime_utr_end:
                    if file_type == "gtf":
                        five_prime_utr_attribute = (
                            f'gene_id "{gene.iloc[0]["gene_id"]}"; '
                            f'Parent "rna-{transcript_id}"; '
                            f'five_prime_utr_id "tss_id-{tss_group.iloc[0]["tss_id"]}"; '
                            f'transcript_id "{transcript_id}"'
                        )
                    elif file_type == "gff":
                        five_prime_utr_attribute = (
                            f'ID=five_prime_utr-tss_id-{tss_group.iloc[0]["tss_id"]};'
                            f'gene_id={gene.iloc[0]["gene_id"]};'
                            f'Parent=rna-{transcript_id};'
                            f'transcript_id={transcript_id}'
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

        # 複数のtranscript_idを持つ場合の処理
        else:
            # 複数のtranscript_idを持つ場合、TSSに最も近いtranscript_idを選択
            print(f"Gene ID {gene_id} has multiple transcript_ids: {transcript_ids}")
            tss_pos = tss_group['tss_pos'].iloc[0]

            # 各transcriptのCDS開始点または終了点を取得し、TSSとの距離を計算
            tss_to_exon_distance = {}  # 修正: transcript_distances → tss_to_cds_distance
            for transcript_id in transcript_ids:
                transcript_exon = exons[exons['transcript_id'] == transcript_id]
                if transcript_exon.empty:
                    continue
                if gene.iloc[0]['strand'] == '+':
                    exon_start = transcript_exon['start'].min()
                else:
                    exon_start = transcript_exon['end'].max()
                tss_to_exon_distance[transcript_id] = abs(exon_start - tss_pos)

            # 最小値を取得
            min_distance = min(tss_to_exon_distance.values())

            # 最小値を持つすべての transcript_id を取得
            closest_transcripts = [tid for tid, distance in tss_to_exon_distance.items() if distance == min_distance]

            # 最小値を持つすべての transcript_id に対して処理を実行
            for transcript_id in closest_transcripts:
                print(f"Processing transcript_id {transcript_id} for gene_id {gene_id}")

                # 選ばれた transcript_id の情報を取得
                transcript_cds = cds[cds['transcript_id'] == transcript_id]

                # 以下、選択された transcript_id を用いて処理を実行
                strand = gene.iloc[0]['strand']

                # 選択された transcript_id の start & end を取得 (mRNA, transcript feature が含まれないgtfに対応するため)
                transcript_start = exons[exons['transcript_id'] == transcript_id]['start'].iloc[0]
                transcript_end = exons[exons['transcript_id'] == transcript_id]['end'].iloc[0]

                # gene の開始または終了位置を更新
                if strand == '+':
                    if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
                        # transcript start が gene start と一致している場合のみ更新
                        if transcript_start == gene.iloc[0]['start']:
                            new_gene_start = int(tss_pos)
                            gene.iloc[0, gene.columns.get_loc('start')] = new_gene_start
                            # 更新された gene 行を updated_gtf に追加
                            updated_gene = gene.iloc[0].to_dict()
                            updated_gene['original_index'] = gene.index[0]
                            updated_gtf.append(updated_gene)
                            print(f"Updated gene start for gene_id {gene_id} to {tss_pos} (matched transcript start).")
                        else:
                            print(f"Gene start for gene_id {gene_id} not updated (transcript start does not match gene start).")
                            original_gene = gene.iloc[0].to_dict()
                            original_gene['original_index'] = gene.index[0]
                            unchanged_gtf.append(original_gene)
                            print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")
                            print(original_gene)
                    else:
                        print(f"Gene start for gene_id {gene_id} not updated (TSS position is inside CDS).")
                        # オリジナルの gene 行を updated_gtf に追加
                        original_gene = gene.iloc[0].to_dict()
                        original_gene['original_index'] = gene.index[0]
                        unchanged_gtf.append(original_gene)
                        print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")

                    # transcript の更新処理
                    if not transcripts.empty:
                        updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                        if not updated_transcript.empty:
                            if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
                                updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(tss_pos)
                                updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                                print(f"Updated transcript for transcript_id {transcript_id}:")
                                print(updated_transcript.iloc[0])
                            else:
                                print(f"Transcript start for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                                # オリジナルの transcript or mRNA 行を updated_gtf に追加
                                updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].iloc[0].to_dict()
                                updated_transcript['original_index'] = transcripts[transcripts['transcript_id'] == transcript_id].index[0]
                                updated_gtf.append(updated_transcript)
                                print(f"Original transcript entry added to updated_gtf for transcript_id {transcript_id}:")
    
                        # 選択されなかった transcript_id の行をそのまま追加
                        other_transcripts = transcripts[transcripts['transcript_id'] != transcript_id]
                        for idx, transcript_row in other_transcripts.iterrows():
                            unchanged_gtf.append({**transcript_row.to_dict(), 'original_index': idx})
                            ###print(f"Unmodified transcript added for transcript_id {transcript_row['transcript_id']}:")
                            ###print(transcript_row)

                    # exon の更新処理
                    related_exons = exons[exons['attribute'].str.contains(
                        rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                    )]
                    if not related_exons.empty:
                        if tss_pos < cds[cds['transcript_id'] == transcript_id]['start'].min():
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
                                ###print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                                ###print(exon_row)
                        else:
                            print(f"Exon start for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                            # オリジナルの exon 行を updated_gtf に追加
                            original_exon = related_exons.iloc[0].to_dict()
                            original_exon['original_index'] = related_exons.index[0]
                            unchanged_gtf.append(original_exon)
                            print(f"Original exon entry added to updated_gtf for transcript_id {transcript_id}:")
                        
                        # 選ばれなかった transcript_id を持つ exon をそのまま追加
                        unselected_exons = exons[~exons['attribute'].str.contains(
                            rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                        )]
                        for idx, exon_row in unselected_exons.iterrows():
                            unchanged_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                            ###print(f"Unmodified exon added for unselected transcript_id {exon_row['attribute']}:")
                            ###print(exon_row)

                    # five_prime_UTR の計算と追加処理    
                    # 最初の CDS を基準に five_prime_utr_end を計算
                    five_prime_utr_start = tss_pos  # 修正: new_gene_start → tss_pos
                    five_prime_utr_end = transcript_cds['start'].min() - 1
            
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

                # strand == '-' の場合の処理
                else:
                    if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
                        # TSSに最も近いtranscriptのendがgene_endと一致する場合のみgene_endを更新
                        if transcript_end == gene.iloc[0]['end']:  
                            new_gene_end = int(tss_pos)  # 更新されたgene_endをnew_gene_endに設定
                            gene.iloc[0, gene.columns.get_loc('end')] = new_gene_end
                            print(f"Updated gene end for gene_id {gene_id} to {tss_pos} (matched transcript end).")
                            # 更新された gene 行を updated_gtf に追加
                            updated_gene = gene.iloc[0].to_dict()
                            updated_gene['original_index'] = gene.index[0]
                            updated_gtf.append(updated_gene)
                        else:
                            print(f"Gene start for gene_id {gene_id} not updated (transcript start does not match gene start).")
                            original_gene = gene.iloc[0].to_dict()
                            original_gene['original_index'] = gene.index[0]
                            unchanged_gtf.append(original_gene)
                            print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")
                            print(original_gene)    
                    else:
                        print(f"Gene start for gene_id {gene_id} not updated (TSS position is inside CDS).")
                        # オリジナルの gene 行を updated_gtf に追加
                        original_gene = gene.iloc[0].to_dict()
                        original_gene['original_index'] = gene.index[0]
                        unchanged_gtf.append(original_gene)
                        print(f"Original gene entry added to updated_gtf for gene_id {gene_id}:")

                    # transcript の更新処理を追加
                    if not transcripts.empty:
                        updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
                        if not updated_transcript.empty:
                            if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
                                if not updated_transcript.empty:
                                    updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(tss_pos)
                                    updated_gtf.append({**updated_transcript.iloc[0].to_dict(), 'original_index': updated_transcript.index[0]})
                                    print(f"Updated transcript for transcript_id {transcript_id}:")
                                    print(updated_transcript.iloc[0])
                            else:
                                print(f"Transcript end for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                                # オリジナルの transcript or mRNA 行を updated_gtf に追加
                                updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].iloc[0].to_dict()
                                updated_transcript['original_index'] = transcripts[transcripts['transcript_id'] == transcript_id].index[0]
                                unchanged_gtf.append(updated_transcript)
                                print(f"Original transcript entry added to updated_gtf for transcript_id {transcript_id}:")
    
                        # 選択されなかった transcript_id の行をそのまま追加
                        other_transcripts = transcripts[transcripts['transcript_id'] != transcript_id]
                        for idx, transcript_row in other_transcripts.iterrows():
                            unchanged_gtf.append({**transcript_row.to_dict(), 'original_index': idx})
                            ###print(f"Unmodified transcript added for transcript_id {transcript_row['transcript_id']}:")
                            ###print(transcript_row)

                    # exon の更新処理を追加
                    # 選択された transcript_id に関連する exon を取得
                    related_exons = exons[exons['attribute'].str.contains(
                        rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                    )]
                    if not related_exons.empty:
                        if tss_pos > cds[cds['transcript_id'] == transcript_id]['end'].max():
                            # 更新された exon を追加
                            largest_end_exon = related_exons.loc[related_exons['end'].idxmax()].copy()
                            largest_end_exon['end'] = int(tss_pos)  # TSS に基づいて end を更新
                            updated_gtf.append({**largest_end_exon.to_dict(), 'original_index': largest_end_exon.name})
                            print(f"Updated exon for transcript_id {transcript_id} (largest end):")
                            print(largest_end_exon)
                            # 他の exon をそのまま追加
                            other_exons = related_exons[related_exons.index != largest_end_exon.name]
                            for idx, exon_row in other_exons.iterrows():
                                updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                                ###print(f"Unmodified exon added for transcript_id {exon_row['attribute']}:")
                                ###print(exon_row)
                        else:
                            print(f"Exon end for transcript_id {transcript_id} not updated (TSS position is inside CDS).")
                            # オリジナルの exon 行を updated_gtf に追加
                            original_exon = related_exons.iloc[0].to_dict()
                            original_exon['original_index'] = related_exons.index[0]
                            updated_gtf.append(original_exon)
                            print(f"Original exon entry added to updated_gtf for transcript_id {transcript_id}:")
                        # 選ばれなかった transcript_id を持つ exon をそのまま追加
                        unselected_exons = exons[~exons['attribute'].str.contains(
                            rf'transcript_id[ =]"?{transcript_id}"?', regex=True
                        )]
                        for idx, exon_row in unselected_exons.iterrows():
                            updated_gtf.append({**exon_row.to_dict(), 'original_index': idx})
                            ###print(f"Unmodified exon added for unselected transcript_id {exon_row['attribute']}:")
                            ###print(exon_row)

                    # five_prime_UTR の計算と追加処理    
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

    # unchanged_gtf から更新された行のオリジナル情報を削除
    try:
        # updated_gtf の original_index を取得
        updated_indices = {entry['original_index'] for entry in updated_gtf if 'original_index' in entry}

        # unchanged_gtf をフィルタリングして、updated_gtf に含まれる original_index を持つ行を削除
        unchanged_gtf = [entry for entry in unchanged_gtf if entry.get('original_index') not in updated_indices]

        print("Removed updated entries from unchanged_gtf.")
    except Exception as e:
        print(f"Error while removing updated entries from unchanged_gtf: {e}")
        raise

    # unchanged_gtf をファイルに出力
    try:
        unchanged_gtf_df = pd.DataFrame(unchanged_gtf)
        if 'original_index' in unchanged_gtf_df.columns:
            unchanged_gtf_df = unchanged_gtf_df.sort_values(by=['original_index'], ascending=True)
        unchanged_gtf_output_file = "unchanged_gtf_debug_output.tsv"
        unchanged_gtf_df.to_csv(unchanged_gtf_output_file, sep='\t', index=False)
        print(f"Unchanged GTF data has been written to {unchanged_gtf_output_file}")
    except Exception as e:
        print(f"Error while writing unchanged_gtf to file: {e}")
        raise

    # updated_gtf をファイルに出力
    try:
        updated_gtf_df = pd.DataFrame(updated_gtf)
        if 'original_index' in updated_gtf_df.columns:
            updated_gtf_df = updated_gtf_df.sort_values(by=['original_index'], ascending=True)
        updated_gtf_output_file = "updated_gtf_debug_output.tsv"
        updated_gtf_df.to_csv(updated_gtf_output_file, sep='\t', index=False)
        print(f"Updated GTF data has been written to {updated_gtf_output_file}")
    except Exception as e:
        print(f"Error while writing updated_gtf to file: {e}")
        raise

    # updated_gtf, unchanged_gtf, cds_gtf を結合
    try:
        final_gtf_df = pd.DataFrame(updated_gtf + unchanged_gtf + cds_gtf)
    except Exception as e:
        print(f"Error while combining updated_gtf, unchanged_gtf, and cds_gtf: {e}")
        raise

    # デバッグ: 結合後のデータフレームを確認
    print("Debug: Combined DataFrame before sorting:")
    print(final_gtf_df.head(10))

    # original_index と feature でソート
    if 'original_index' in final_gtf_df.columns and 'feature' in final_gtf_df.columns:
        try:
            # feature列のカスタム順序を設定
            feature_order = ['region', 'gene', 'pseudogene', 'transcript', 'mRNA', 'lnc_RNA', 'tRNA', 'miRNA', 'primary_transcript', 'rRNA', 
                            'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon', 
                            'V_gene_segment', 'C_gene_segment', 'RNase_P_RNA', 'antisense_RNA', 'match', 'cDNA_match']
            final_gtf_df['feature'] = pd.Categorical(final_gtf_df['feature'], categories=feature_order, ordered=True)

            # original_index と feature でソート
            final_gtf_df = final_gtf_df.sort_values(by=['original_index', 'feature'], ascending=[True, True])
        except Exception as e:
            print(f"Error while sorting by 'original_index' and 'feature': {e}")
            raise
    else:
        print("Warning: 'original_index' or 'feature' column is missing. Skipping sorting by these columns.")

    # 重複行を削除
    try:
        final_gtf_df = final_gtf_df.drop_duplicates()
    except Exception as e:
        print(f"Error while dropping duplicates: {e}")
        raise

    # デバッグ: 重複削除後のデータフレームを確認
    print("Debug: Final DataFrame after deduplication:")
    print(final_gtf_df.head(10))

    # gtf_df と final_gtf_df の feature 列の数を比較してサマライズ
    try:
        # gtf_df の feature 別のカウント
        gtf_feature_counts = gtf_df_original['feature'].value_counts()

        # final_gtf_df の feature 別のカウント
        final_feature_counts = final_gtf_df['feature'].value_counts()

        # デバッグメッセージとして出力
        print("Feature counts in gtf_df:")
        print(gtf_feature_counts)
        print("\nFeature counts in final_gtf_df:")
        print(final_feature_counts)

        # 差分を計算して出力
        try:
            feature_diff = gtf_feature_counts.subtract(final_feature_counts, fill_value=0)
            # 差分を整数型に変換
            feature_diff = feature_diff.astype(int)
            print("\nDifference in feature counts (gtf_df - final_gtf_df):")
            print(feature_diff)

            # final_gtf_df に含まれていない feature 行を抽出
            missing_features = ['gene', 'mRNA', 'transcript']
            if 'original_index' in gtf_df_original.columns and 'original_index' in final_gtf_df.columns:
                missing_rows = gtf_df_original[
                    (gtf_df_original['feature'].isin(missing_features)) &
                    (~gtf_df_original['original_index'].isin(final_gtf_df['original_index']))
                ]
                # デバッグ: 含まれていない行を出力
                print("\nRows in gtf_df_original with features 'gene', 'mRNA', 'transcript' not in final_gtf_df:")
                print(missing_rows)
            else:
                print("Error: 'original_index' column is missing in one of the DataFrames.")
                missing_rows = pd.DataFrame()  # 空のデータフレームを作成

            # デバッグ: 含まれていない行を出力
            print("\nRows in gtf_df_original with features 'gene', 'mRNA', 'transcript' not in final_gtf_df:")
            print(missing_rows)

            # 必要に応じてファイルに出力
            missing_rows_output_file = "missing_features_debug_output.tsv"
            missing_rows.to_csv(missing_rows_output_file, sep='\t', index=False)
            print(f"Missing feature rows have been written to {missing_rows_output_file}")
        except Exception as e:
            print(f"Error while summarizing feature counts or extracting missing rows: {e}")
            raise
    except Exception as e:
        print(f"Error while summarizing feature counts: {e}")
        raise

    # データ部分を文字列形式に変換
    try:
        # 最初の9列のみを保持
        final_gtf_df = final_gtf_df.iloc[:, :9]

        # データを文字列形式に変換
        gtf_data_lines = final_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()
    except Exception as e:
        print(f"Error while converting DataFrame to CSV lines: {e}")
        raise

    # コメント行を先頭に追加し、データ部分を結合
    try:
        output_lines = comment_lines + gtf_data_lines
    except Exception as e:
        print(f"Error while combining comment lines with GTF data: {e}")
        raise

    # ファイルに書き出し
    try:
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_lines) + '\n')
        print("File saved successfully.")
    except Exception as e:
        print(f"Error while writing to output file: {e}")
        raise

    # デバッグ: 最終的なデータフレームのプレビュー
    print("Debug: Final GTF DataFrame preview:")
    print(final_gtf_df.head(50))  # 最初の50行を確認

if __name__ == "__main__":
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description="GTF/GFFファイルをTSS情報で更新するスクリプト")
    parser.add_argument("gtf_file", help="入力GTFまたはGFFファイルのパス")  # 修正
    parser.add_argument("tss_file", help="入力TSSファイルのパス")  # 修正
    parser.add_argument("output_file", help="出力GTFファイルのパス")  # 修正
    args = parser.parse_args()

    # 関数を実行
    update_gtf_with_tss(args.gtf_file, args.tss_file, args.output_file)