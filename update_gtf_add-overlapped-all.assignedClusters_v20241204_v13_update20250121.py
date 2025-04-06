import pandas as pd
import csv

def update_gtf_with_tss(gtf_file, tss_file, output_file):
    # GTFファイルとTSSファイルの読み込み
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    tss_df = pd.read_csv(tss_file, sep='\t', header=0)

    # gene_idとtranscript_idを抽出
    gtf_df['gene_id'] = gtf_df['attribute'].str.extract(r'gene_id "([^"]+)"')
    gtf_df['transcript_id'] = gtf_df['attribute'].str.extract(r'transcript_id "([^"]+)"')
    ## NaNをデフォルト値に置き換え
    #gtf_df['gene_id'] = gtf_df['gene_id'].fillna('unknown_gene')
    #gtf_df['transcript_id'] = gtf_df['transcript_id'].fillna('unknown_transcript')

    tss_gene_ids = tss_df['gene_id'].unique()
    updated_gtf = []

    # TSSのgene_idに基づいて更新
    for gene_id in tss_gene_ids:
        gene_group = gtf_df[gtf_df['gene_id'] == gene_id]
        if gene_group.empty:
            continue

        gene = gene_group[gene_group['feature'] == 'gene']
        exons = gene_group[gene_group['feature'] == 'exon']
        start_codons = gene_group[gene_group['feature'] == 'start_codon']
        stop_codons = gene_group[gene_group['feature'] == 'stop_codon']
        transcripts = gene_group[gene_group['feature'] == 'transcript']
        other_features = gene_group[~gene_group['feature'].isin(['gene', 'exon', 'transcript'])]

        tss_group = tss_df[tss_df['gene_id'] == gene_id]

        if gene.empty or start_codons.empty or stop_codons.empty:
            continue

        for _, start_codon in start_codons.iterrows():
            transcript_id = start_codon['transcript_id']
            strand = gene.iloc[0]['strand']
            
            if strand == '+':
                new_gene_start = tss_group['tss_pos'].iloc[0]
            else:
                new_gene_end = tss_group['tss_pos'].iloc[0]

            if strand == '+':
                exons_for_transcript = exons[exons['transcript_id'] == transcript_id]
                first_exon = exons_for_transcript.sort_values(by='start').iloc[0]
                other_exons = exons_for_transcript[(exons_for_transcript['start'] != first_exon['start']) | (exons_for_transcript['end'] != first_exon['end'])]
                first_exon_start = new_gene_start
                five_prime_utr_start = new_gene_start
                five_prime_utr_end = start_codon['start'] - 1
                # first_exonの開始位置を更新
                first_exon['start'] = first_exon_start
            else:
                exons_for_transcript = exons[exons['transcript_id'] == transcript_id]
                first_exon = exons_for_transcript.sort_values(by='end', ascending=False).iloc[0]
                other_exons = exons_for_transcript[(exons_for_transcript['start'] != first_exon['start']) | (exons_for_transcript['end'] != first_exon['end'])]
                first_exon_end = new_gene_end
                five_prime_utr_start = start_codon['end'] + 1
                five_prime_utr_end = new_gene_end
                # first_exonの終了位置を更新
                first_exon['end'] = first_exon_end

            updated_gene = gene.copy()
            if strand == '+':
                updated_gene.iloc[0, updated_gene.columns.get_loc('start')] = int(new_gene_start)
            else:
                updated_gene.iloc[0, updated_gene.columns.get_loc('end')] = int(new_gene_end)
            updated_gtf.append(updated_gene.to_dict(orient='records')[0])

            # 更新されたfirst_exonを追加
            updated_gtf.append(first_exon.to_dict())

            if not other_exons.empty:
                updated_gtf.extend(other_exons.to_dict(orient='records'))

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
                    'attribute': [gene.iloc[0]['attribute']],
                    'gene_id': [gene_id],  # 追加
                    'transcript_id': [transcript_id]  # 追加
                })
                updated_gtf.append(five_prime_utr.to_dict(orient='records')[0])
            
            # 対応するトランスクリプトの更新
            updated_transcript = transcripts[transcripts['transcript_id'] == transcript_id].copy()
            if strand == '+':
                updated_transcript.iloc[0, updated_transcript.columns.get_loc('start')] = int(new_gene_start)
            else:
                updated_transcript.iloc[0, updated_transcript.columns.get_loc('end')] = int(new_gene_end)
            updated_gtf.append(updated_transcript.to_dict(orient='records')[0])
        
        # 他の要素をそのまま追加
        updated_gtf.extend(other_features.to_dict(orient='records'))

    # デバッグ用：updated_gtfを出力
    print("Final updated_gtf:")
    for row in updated_gtf:
        if 'gene_id' not in row or 'transcript_id' not in row:
            print("Error in row:", row)

    # 未更新の行を追加
    updated_row_ids = set((row['gene_id'], row['transcript_id']) for row in updated_gtf)
    non_updated_rows = gtf_df[
        ~gtf_df.apply(lambda row: (row['gene_id'], row['transcript_id']) in updated_row_ids, axis=1)
    ]
    updated_gtf.extend(non_updated_rows.to_dict(orient='records'))

    # 保存
    updated_gtf_df = pd.DataFrame(updated_gtf)
    updated_gtf_df['attribute'] = updated_gtf_df['attribute'].str.split(';').str[0] + ';'  # 余分な部分を削除
    updated_gtf_df.to_csv(output_file, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

# 関数を実行
update_gtf_with_tss(
    '/Volumes/project＿HtriTSS＿BU/Project_HtriGenePrediction/BRAKER3_rerun_vRMrun202404/braker/braker_correctID_v3_add_gene_id_inFeatureTranscript_by-gemini_20241203.gtf',
    '/Volumes/project＿HtriTSS＿BU/Project_HtriTSS/CAGE_analysis_v20241025_RM-adapterReads_selectedGenomeSeq_gtf-v03/wd4selected-each3gonad_and_tissueMix/all_tss_feature_uniq.gene.tsv',
    'braker_correctID_v3_add_gene_id_inFeatureTranscript_by-gemini_20241204_w_sort.all_tss_feature_uniq.gene_script-v13.gtf'
)