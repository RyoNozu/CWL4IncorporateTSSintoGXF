import polars as pl
import csv

def update_gtf_with_tss(gtf_file, tss_file, output_file):
    # GTFファイルとTSSファイルの読み込み
    gtf_df = pl.read_csv(gtf_file, separator='\t', comment_prefix='#', has_header=False, 
                         new_columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    tss_df = pl.read_csv(tss_file, separator='\t', has_header=True)

    # gene_idとtranscript_idを抽出
    gtf_df = gtf_df.with_columns([
        gtf_df['attribute'].str.extract(r'gene_id "([^"]+)"').alias('gene_id'),
        gtf_df['attribute'].str.extract(r'transcript_id "([^"]+)"').alias('transcript_id')
    ])

    tss_gene_ids = tss_df['gene_id'].unique()
    updated_gtf = []
    updated_row_ids = set()  # ここで初期化

    # TSSのgene_idに基づいて更新
    for gene_id in tss_gene_ids:
        gene_group = gtf_df.filter(gtf_df['gene_id'] == gene_id)
        if gene_group.is_empty():
            continue

        gene = gene_group.filter(gene_group['feature'] == 'gene')
        exons = gene_group.filter(gene_group['feature'] == 'exon')
        start_codons = gene_group.filter(gene_group['feature'] == 'start_codon')
        stop_codons = gene_group.filter(gene_group['feature'] == 'stop_codon')
        transcripts = gene_group.filter(gene_group['feature'] == 'transcript')
        other_features = gene_group.filter(~gene_group['feature'].is_in(['gene', 'exon', 'transcript']))

        tss_group = tss_df.filter(tss_df['gene_id'] == gene_id)

        if gene.is_empty() or start_codons.is_empty() or stop_codons.is_empty():
            continue

        for start_codon in start_codons.iter_rows(named=True):
            transcript_id = start_codon['transcript_id']
            strand = gene[0, 'strand']
            
            if strand == '+':
                new_gene_start = tss_group[0, 'tss_pos']
            else:
                new_gene_end = tss_group[0, 'tss_pos']

            if strand == '+':
                exons_for_transcript = exons.filter(exons['transcript_id'] == transcript_id)
                first_exon = exons_for_transcript.sort('start').row(0)
                other_exons = exons_for_transcript.filter(
                    (exons_for_transcript['start'] != first_exon['start']) | 
                    (exons_for_transcript['end'] != first_exon['end'])
                )
                first_exon_start = new_gene_start
                five_prime_utr_start = new_gene_start
                five_prime_utr_end = start_codon['start'] - 1
                # first_exonの開始位置を更新
                first_exon['start'] = first_exon_start
            else:
                exons_for_transcript = exons.filter(exons['transcript_id'] == transcript_id)
                first_exon = exons_for_transcript.sort('end', reverse=True).row(0)
                other_exons = exons_for_transcript.filter(
                    (exons_for_transcript['start'] != first_exon['start']) | 
                    (exons_for_transcript['end'] != first_exon['end'])
                )
                first_exon_end = new_gene_end
                five_prime_utr_start = start_codon['end'] + 1
                five_prime_utr_end = new_gene_end
                # first_exonの終了位置を更新
                first_exon['end'] = first_exon_end

            updated_gene = gene.clone()
            if strand == '+':
                updated_gene[0, 'start'] = int(new_gene_start)
            else:
                updated_gene[0, 'end'] = int(new_gene_end)
            updated_gtf.append(updated_gene.to_dicts()[0])

            # 更新されたfirst_exonを追加
            updated_gtf.append(first_exon)

            if not other_exons.is_empty():
                updated_gtf.extend(other_exons.to_dicts())

            if five_prime_utr_start <= five_prime_utr_end:
                five_prime_utr = pl.DataFrame({
                    'seqname': [gene[0, 'seqname']],
                    'source': ['.'],
                    'feature': ['five_prime_UTR'],
                    'start': [int(five_prime_utr_start)],
                    'end': [int(five_prime_utr_end)],
                    'score': ['.'],
                    'strand': [strand],
                    'frame': ['.'],
                    'attribute': [gene[0, 'attribute']],
                    'gene_id': [gene_id],  # 追加
                    'transcript_id': [transcript_id]  # 追加
                })
                updated_gtf.append(five_prime_utr.to_dicts()[0])
            
            # 対応するトランスクリプトの更新
            updated_transcript = transcripts.filter(transcripts['transcript_id'] == transcript_id).clone()
            if strand == '+':
                updated_transcript[0, 'start'] = int(new_gene_start)
            else:
                updated_transcript[0, 'end'] = int(new_gene_end)
            updated_gtf.append(updated_transcript.to_dicts()[0])

            # updated_row_ids に追加
            updated_row_ids.add((gene_id, transcript_id))
        
        # 他の要素をそのまま追加
        updated_gtf.extend(other_features.to_dicts())

    # デバッグ用：updated_gtfを出力
    print("Final updated_gtf:")
    for row in updated_gtf:
        if 'gene_id' not in row or 'transcript_id' not in row:
            print("Error in row:", row)

    # updated_row_ids を DataFrame に変換（データ型を明示的に指定）
    updated_row_ids_df = pl.DataFrame(
        list(updated_row_ids),
        schema={'gene_id': pl.Utf8, 'transcript_id': pl.Utf8}
    )

    print(gtf_df.dtypes)
    print(updated_row_ids_df.dtypes)

    # gtf_df のデータ型を統一
    gtf_df = gtf_df.with_columns([
        gtf_df['gene_id'].cast(pl.Utf8),
        gtf_df['transcript_id'].cast(pl.Utf8)
    ])

    # gtf_df に対してフィルタリング
    non_updated_rows = gtf_df.join(updated_row_ids_df, on=['gene_id', 'transcript_id'], how='anti')

    # 未更新の行を追加
    updated_gtf.extend(non_updated_rows.to_dicts())

    # 保存
    updated_gtf_df = pl.DataFrame(updated_gtf)

    # attribute 列の処理
    updated_gtf_df = updated_gtf_df.with_columns([
        updated_gtf_df['attribute'].map_elements(
            lambda x: x.split(';')[0] + ';' if isinstance(x, str) else x,
            return_dtype=pl.Utf8
        ).alias('attribute')
    ])

    print(updated_gtf_df['attribute'].head())
    print(updated_gtf_df['attribute'].dtype)

    # ファイルに書き出し
    updated_gtf_df.write_csv(output_file, separator='\t', include_header=False)

# 関数を実行
update_gtf_with_tss(
    '/Volumes/project＿HtriTSS＿BU/genomes_NCBI/GCF_000002035.6_GRCz11_acc20230614/4script-check_filtered-NW.GCF_000002035.6_GRCz11_genomic.gtf',
    '/Volumes/project＿HtriTSS＿BU/Project_HtriTSS/CAGE_analysis_v20250108_RM-adapterReads_zebrafish-PRJNA799647/all_tss_feature_uniq.gene.tsv',
    'update_by-polars.gtf'
)