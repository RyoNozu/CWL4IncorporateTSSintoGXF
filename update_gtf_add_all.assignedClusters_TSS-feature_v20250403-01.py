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
        cds = gene_group[gene_group['feature'] == 'CDS']
        other_features = gene_group[~gene_group['feature'].isin(['gene', 'exon', 'transcript', 'CDS'])]

        tss_group = tss_df[tss_df['gene_id'] == gene_id]

        if gene.empty or cds.empty:
            print(f"Skipping gene_id {gene_id} due to missing gene or CDS.")
            continue

        for _, start_codon_row in cds.iterrows():
            transcript_id = start_codon_row['transcript_id']
            if pd.isna(transcript_id):
                print(f"Transcript ID is missing for gene_id {gene_id}. Skipping.")
                continue

            print(f"Processing transcript_id: {transcript_id}")
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
            for transcript_id in cds['transcript_id'].unique():
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
                    five_prime_utr = pd.DataFrame({
                        'seqname': [gene.iloc[0]['seqname']],
                        'source': ['.'],
                        'feature': ['five_prime_UTR'],
                        'start': [int(five_prime_utr_start)],
                        'end': [int(five_prime_utr_end)],
                        'score': ['.'],
                        'strand': [strand],
                        'frame': ['.'],
                        'attribute': [f'transcript_id "{transcript_id}"; ' + gene.iloc[0]['attribute']]
                    })
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

    # 最終的に DataFrame に変換
    updated_gtf_df = pd.DataFrame(updated_gtf)

    # 元の順番を保持するためにソート
    updated_gtf_df = updated_gtf_df.sort_values(by='original_index').drop(columns=['original_index'])

    # 10列目以降を削除
    updated_gtf_df = updated_gtf_df.iloc[:, :9]

    # 重複を取り除く
    updated_gtf_df = updated_gtf_df.drop_duplicates()

    # コメント行を元の位置に挿入
    output_lines = lines.copy()  # 元の行をコピー

    # コメント行以外のデータを文字列形式に変換して追加
    gtf_data_lines = updated_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()
    for i, line in enumerate(output_lines):
        if not line.startswith('#'):
            # コメント行以外の部分を置き換え
            output_lines[i:] = gtf_data_lines
            break

    # 保存処理
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')

if __name__ == "__main__":
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description="GTF/GFFファイルをTSS情報で更新するスクリプト")
    parser.add_argument("gtf_file", help="入力GTFまたはGFFファイルのパス")
    parser.add_argument("tss_file", help="入力TSSファイルのパス")
    parser.add_argument("output_file", help="出力GTFファイルのパス")
    args = parser.parse_args()

    # 関数を実行
    update_gtf_with_tss(args.gtf_file, args.tss_file, args.output_file)