#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "step5: The process of updating the GFF format file from identifying TSS (transcription start sites) from CAGE-seq data"
doc: |
  "
  The process of updating the GFF format file from identifying TSS - transcription start sites - from paired-end CAGE-seq data.
  This workflow consists of the following files:
  (1) Tools/06_combined_exec_TSSr.cwl,
  (2) Tools/07_join_all_assignedClusters.cwl,
  (3) Tools/08_uniq_tss_feature.cwl,
  (4) Tools/09_update_gtf.cwl
  "

inputs:
  # inputs parameters for Tools/06_combined_exec_TSSr.cwl
  - id: reference_gtf_file
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "GFF format reference genome annotation file"
    doc: "GFF format reference genome annotation file"

  - id: seed_file
    type: File
    format: edam:format_3671 # text
    label: "seed file"
    doc: "seed file for BSgenome package"

  - id: threads
    type: int
    label: "number of threads"
    doc: "number of threads"
    default: 16

  - id: input_file_type
    type: string
    label: "input files type"
    doc: "input BAM files type (bamPairedEnd or bam). For example, if the processed CAGE-seq files are paired-end files, select “bamPairedEnd”."
    default: "bamPairedEnd"

  - id: organism_name
    type: string
    label: "organism name"
    doc: "organism name"
    default: "Halichoeres trimaculatus"

  - id: metadata_file
    type: File
    format: edam:format_3752 # csv
    label: "metadata file"
    doc: "metadata file for grouping bam files"

  - id: genome_seqs_dir
    type: Directory
    label: "genome sequences directory"
    doc: "genome sequences directory processed by seqkit process"

  - id: annotation_region_upstream
    type: int
    label: "annotation region upstream"
    doc: "The maximum allowable distance for assigning a dominant TSS (transcription start site) to its downstream gene, when the TSS is located upstream of 5'-end of the CDS or the the transcript. Default: 1000 bp"
    default: 1000

  - id: annotation_region_downstream
    type: int
    label: "annotation region downstream"
    doc: "The maximum allowable distance for assigning a TSS when the dominant TSS is located downstream of 5'-end of the CDS or the the transcript. Default: 0 (i.e., TSSs located downstream of the 5'-end of feature (CDS or transcript) are not assigned to the gene)"
    default: 0

  - id: annotation_type
    type: string
    label: "annotation type"
    doc: "annotation type (\"genes\" or \"transcripts\")"
    default: "transcripts"

  - id: bam_files
    type: File[]
    format: edam:format_2572 # bam file
    label: "bam files"
    doc: "bam files processed by STAR mapping process"

  # inputs parameters for Tools/09_update_gtf.cwl
  - id: update_gtf_filename
    type: string
    label: "File name in GFF format"
    doc: "File name in GFF format with TSS (transcription start sites) information added"
    default: "add_tss_feature.gtf"

steps:
  # Tools/06_combined_exec_TSSr.cwl
  - id: exec_TSSr
    run: ../Tools/06_combined_exec_TSSr.cwl
    in:
      referenceFile: reference_gtf_file
      seedFile: seed_file
      threads: threads
      inputFilesType: input_file_type
      organismName: organism_name
      metadata: metadata_file
      genome_seqs_dir: genome_seqs_dir
      annotation_region_upstream: annotation_region_upstream
      annotation_region_downstream: annotation_region_downstream
      annotation_type: annotation_type
      bam_files: bam_files
    out:
      - clustered_consensus_TSSs_bed
      - clustered_consensus_TSSs_txt
      - clustered_assigned_TSSs_txt

  # Tools/07_join_all_assignedClusters.cwl
  - id: join_all_assignedClusters
    run: ../Tools/07_join_all_assignedClusters.cwl
    in:
      assigned_clusters_files: exec_TSSr/clustered_assigned_TSSs_txt
    out:
      - joined_clusters

  # Tools/08_uniq_tss_feature.cwl
  - id: uniq_tss_feature
    run: ../Tools/08_uniq_tss_feature.cwl
    in:
      all_joined_assigned_clusters_file: join_all_assignedClusters/joined_clusters
    out:
      - all_cage_cluster_feature_uniq_gene_file
      - all_tss_feature_uniq_gene_file
      - all_tss_feature_file

  # Tools/09_update_gtf.cwl
  - id: update_gtf
    run: ../Tools/09_update_gtf.cwl
    in:
      gtf_file: reference_gtf_file
      tss_file: uniq_tss_feature/all_tss_feature_uniq_gene_file
      update_gtf_filename: update_gtf_filename
    out:
      - output_gtf_file
      - stdout_log

outputs:
  # outputs parameters for Tools/06_combined_exec_TSSr.cwl
  - id: clustered_consensus_TSSs_bed
    type: File[]
    format: edam:format_2572 # bed file
    label: "clustered consensus TSSs bed file"
    doc: "clustered consensus TSSs bed file"
    outputSource: exec_TSSr/clustered_consensus_TSSs_bed

  - id: clustered_consensus_TSSs_txt
    type: File[]
    format: edam:format_3671 # text
    label: "clustered consensus TSSs txt file"
    doc: "clustered consensus TSSs txt file"
    outputSource: exec_TSSr/clustered_consensus_TSSs_txt

  - id: clustered_assigned_TSSs_txt
    type: File[]
    format: edam:format_3671 # text
    label: "clustered assigned TSSs txt file"
    doc: "clustered assigned TSSs txt file"
    outputSource: exec_TSSr/clustered_assigned_TSSs_txt
 
  # outputs parameters for Tools/07_join_all_assignedClusters.cwl
  - id: joined_clusters
    type: File
    format: edam:format_3475 # tsv
    label: "joined assigned clusters"
    doc: "joined assigned clusters file containing merged data from all input files"
    outputSource: join_all_assignedClusters/joined_clusters

  # outputs parameters for Tools/08_uniq_tss_feature.cwl
  - id: all_cage_cluster_feature_uniq_gene_file
    type: File
    format: edam:format_3475 # tsv
    label: "all cage cluster feature uniq gene file"
    doc: "Contains unique cluster information for each gene, extracted from all_tss_feature.tsv."
    outputSource: uniq_tss_feature/all_cage_cluster_feature_uniq_gene_file

  - id: all_tss_feature_uniq_gene_file
    type: File
    format: edam:format_3475 # tsv
    label: "all tss feature uniq gene file"
    doc: "Contains unique TSS information for each gene, extracted from all_tss_feature.tsv"
    outputSource: uniq_tss_feature/all_tss_feature_uniq_gene_file

  - id: all_tss_feature_file
    type: File
    format: edam:format_3475 # tsv
    label: "all tss feature uniq gene file"
    doc: "Contains TSS cluster information filtered based on the tag accumulation at the dominant TSS (tags.dominant_tss). Includes clusters from groups with the highest tag accumulation."
    outputSource: uniq_tss_feature/all_tss_feature_file

  # outputs parameters for Tools/09_update_gtf.cwl
  - id: output_gtf_file
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "updated GFF format file"
    doc: "updated GFF format file with TSS (transcription start sites) information added"
    outputSource: update_gtf/output_gtf_file

  - id: stdout_log
    type: File
    label: "stdout log"
    doc: "stdout log"
    outputSource: update_gtf/stdout_log


$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/