import os
import pandas as pd
import argparse
import csv

def extract_gene_id(attribute_column, file_type):
    if file_type == "gtf":
        return attribute_column.str.extract(r'gene_id "([^"]+)"')
    elif file_type == "gff":
        return attribute_column.str.extract(r'(?:ID=gene-|Parent=gene-|gene=)([^;]+)')
    else:
        raise ValueError("Unsupported file type. Only GTF and GFF are supported.")

def extract_transcript_id(attribute_column, file_type):
    if file_type == "gtf":
        return attribute_column.str.extract(r'transcript_id "([^"]+)"')
    elif file_type == "gff":
        return attribute_column.str.extract(r'(?:ID=rna-|Parent=rna-|transcript_id=)([^;]+)')
    else:
        raise ValueError("Unsupported file type. Only GTF and GFF are supported.")

def update_gene_with_tss_and_utr(gene_group, tss_group):
    updated_rows = []
    utr_rows = []

    for _, tss_row in tss_group.iterrows():
        tss_start = tss_row['start']
        tss_end = tss_row['end']

        for _, row in gene_group.iterrows():
            updated_row = row.copy()

            if row['feature'] == 'gene':
                updated_row['start'] = tss_start
            elif row['feature'] in ['transcript', 'mRNA']:
                updated_row['start'] = tss_start
            elif row['feature'] == 'CDS':
                if row['strand'] == '+':
                    utr_start = tss_start
                    utr_end = row['start'] - 1
                else:
                    utr_start = row['end'] + 1
                    utr_end = tss_start

                if utr_start < utr_end:
                    utr_row = row.copy()
                    utr_row['feature'] = "5'UTR"
                    utr_row['start'] = utr_start
                    utr_row['end'] = utr_end
                    utr_rows.append(utr_row)

            updated_rows.append(updated_row)

    return updated_rows + utr_rows

def update_gtf_with_tss(gtf_file, tss_file, output_file):
    file_extension = os.path.splitext(gtf_file)[1].lower()
    if file_extension == ".gtf":
        file_type = "gtf"
    elif file_extension in [".gff", ".gff3"]:
        file_type = "gff"
    else:
        raise ValueError("Unsupported file format. Please provide a .gtf or .gff file.")

    with open(gtf_file, 'r') as f:
        lines = f.readlines()

    comment_lines = [(i, line) for i, line in enumerate(lines) if line.startswith('#')]

    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    tss_df = pd.read_csv(tss_file, sep='\t', header=0)

    gtf_df['gene_id'] = extract_gene_id(gtf_df['attribute'], file_type)
    gtf_df['transcript_id'] = extract_transcript_id(gtf_df['attribute'], file_type)

    gtf_df['original_index'] = gtf_df.index

    updated_gtf = []

    tss_gene_ids = tss_df['gene_id'].unique()
    unprocessed_rows = gtf_df[~gtf_df['gene_id'].isin(tss_gene_ids)]

    for gene_id, gene_group in gtf_df.groupby('gene_id'):
        if gene_id not in tss_gene_ids:
            continue

        transcripts = gene_group[gene_group['feature'].isin(['transcript', 'mRNA'])]
        if transcripts.empty:
            updated_gtf.extend(gene_group.to_dict('records'))
            continue

        if len(transcripts['transcript_id'].unique()) == 1:
            tss_group = tss_df[tss_df['gene_id'] == gene_id]
            if not tss_group.empty:
                updated_gtf.extend(update_gene_with_tss_and_utr(gene_group, tss_group))
            else:
                updated_gtf.extend(gene_group.to_dict('records'))
            continue

        original_tss = gene_group[gene_group['feature'] == 'gene']['start'].min()
        closest_transcript = None
        closest_distance = float('inf')

        for transcript_id, transcript_group in transcripts.groupby('transcript_id'):
            transcript_tss = transcript_group['start'].min()
            distance = abs(transcript_tss - original_tss)
            if distance < closest_distance:
                closest_distance = distance
                closest_transcript = transcript_id

        selected_transcript_group = gene_group[gene_group['transcript_id'] == closest_transcript]

        original_gene_length = gene_group['end'].max() - gene_group['start'].min()
        updated_gene_length = selected_transcript_group['end'].max() - selected_transcript_group['start'].min()

        if updated_gene_length > original_gene_length:
            tss_group = tss_df[tss_df['gene_id'] == gene_id]
            if not tss_group.empty:
                updated_gtf.extend(update_gene_with_tss_and_utr(selected_transcript_group, tss_group))
            else:
                updated_gtf.extend(gene_group.to_dict('records'))
        else:
            updated_gtf.extend(gene_group.to_dict('records'))

    for idx, row in unprocessed_rows.iterrows():
        updated_gtf.append({**row.to_dict(), 'original_index': idx})

    updated_gtf_df = pd.DataFrame(updated_gtf)
    updated_gtf_df = updated_gtf_df.sort_values(by='original_index').drop(columns=['original_index'])
    updated_gtf_df = updated_gtf_df.iloc[:, :9]
    updated_gtf_df = updated_gtf_df.drop_duplicates()

    output_lines = lines.copy()
    gtf_data_lines = updated_gtf_df.to_csv(sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE).splitlines()
    for i, line in enumerate(output_lines):
        if not line.startswith('#'):
            output_lines[i:] = gtf_data_lines
            break

    output_lines = [line.strip() for line in output_lines if line.strip()]

    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GTF/GFFファイルをTSS情報で更新し、5'UTRを追加するスクリプト")
    parser.add_argument("gtf_file", help="入力GTFまたはGFFファイルのパス")
    parser.add_argument("tss_file", help="入力TSSファイルのパス")
    parser.add_argument("output_file", help="出力GTFファイルのパス")
    args = parser.parse_args()

    update_gtf_with_tss(args.gtf_file, args.tss_file, args.output_file)