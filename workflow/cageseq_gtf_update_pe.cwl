#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "CWL4IncorporateTSSintoGXF"
doc: "CWL4IncorporateTSSintoGXF: mapping CAGE-seq data and add TSS information to GXF format (GFF/GTF) file"

requirements:
  # https://www.commonwl.org/user_guide/topics/workflows.html#scattering-steps
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}


inputs:
  #### workflow common inputs (upper-case) ####
  - id: THREADS
    type: int
    label: "threads"
    doc: "threads"
    default: 16

  - id: GENOME_FILE
    type: File
    format: edam:format_3989 # Gzipped format
    label: "genome sequence file (gzipped fasta format) for mapping"
    doc: "genome sequence file (gzipped fasta format) for mapping"


  - id: REFERENCE_GENOME_ANNOTATION
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "GFF format reference genome annotation file"
    doc: "GFF format reference genome annotation file"



  #### fastp unique inputs (paired-end) (lower-case) ####
  - id: fastq1_files
    type: File[]
    format: edam:format_1930
    label: "fastq1 files"
    doc: "fastq1 files (e.g. MK.F1_R1.fastq.gz)"

  - id: fastq2_files
    type: File[]
    format: edam:format_1930
    label: "fastq2 files"
    doc: "fastq2 files (e.g. MK.F1_R2.fastq.gz)"


#### TSSr unique inputs (lower-case) ####
  - id: seed_file
    type: File
    format: edam:format_3671 # text
    label: "seed file"
    doc: "seed file for BSgenome package"

  - id: organism_name
    type: string
    label: "organism name"
    doc: "organism name"

  - id: metadata_file
    type: File
    format: edam:format_3752 # csv
    label: "metadata file"
    doc: "metadata file for grouping bam files"

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
    doc: "annotation type ('genes' or 'transcripts')"
    default: "genes"
    
  - id: update_gtf_filename
    type: string
    label: "File name in GFF format"
    doc: "File name in GFF format with TSS (transcription start sites) information added"
    default: "add_tss_feature.gtf"


steps:
  #### 1. split genome sequences by scaffold ####
  - id: split_genome_seqs
    run: ../Tools/01_split_genome_seqs.cwl
    label: "step1: split genome seqs by scaffold"
    doc: "split genome seqs by scaffold. using seqkit version 2.8.2."
    in:
      genome_seqs: GENOME_FILE
      # output_dir_name: variable
    out: 
      - output_dir
      - split_seqs


  #### 2. trimming fastq files ####
  - id: trimming_fastp_pe
    run: ./01_trimming_fastq_subworkflow_pe.cwl
    label: "step2: trimming fastq files (paired-end)"
    doc: "multiple fastq files trimming process using fastp version 0.23.4 and scatter feature requirement"
    in:
      fastq1_files: fastq1_files
      fastq2_files: fastq2_files
      # compression_level: fastp_compression_level
      threads: THREADS
    out: 
      - trimmed_fastq1_files
      - trimmed_fastq2_files
      - trimming_report_html_files
      - trimming_report_json_files
      - stdout_log_files
      - stderr_log_files


  #### 3. create STAR index ####
  - id: create_star_index
    run: ./02_star_index_subworkflow.cwl
    label: "step3: create STAR index"
    doc: "create STAR index for mapping CAGE-Seq data (step 1: decompress reference genome fasta file, step 2: create STAR index)"
    in:
      create_star_index_threads: THREADS
      input_compressed_file: GENOME_FILE
      # star_index_dir_name: star_index_dir_name
      reference_genome_annotation: REFERENCE_GENOME_ANNOTATION
    out: 
      - star_index_dir
      - stdout_log
      - stderr_log


  #### 4. mapping CAGEseq data####
  - id: mapping_cageseq_data
    run: ./03_star4cageseq_analysis_subworkflow_pe.cwl
    label: "step4: STAR for CAGE-Seq reads analysis"
    doc: "STAR for CAGE-Seq reads analysis"
    in:
      star_index_dir: create_star_index/star_index_dir
      cage_seq_reads_1: trimming_fastp_pe/trimmed_fastq1_files
      cage_seq_reads_2: trimming_fastp_pe/trimmed_fastq2_files
      star_threads: THREADS
    out: 
      - aligned_bam_files
      - final_log_files
      - main_log_files
      - progress_log_files
      - sj_tab_files
      - log_stdout_files
      - log_stderr_files


  #### 5. update GTF file ####
  - id: update_gtf
    run: ./04_tssr_subworkflow_pe.cwl
    label: "step5: The process of updating the GFF format file from identifying TSS (transcription start sites) from CAGE-seq data"
    doc: "The process of updating the GFF format file from identifying TSS - transcription start sites - from paired-end CAGE-seq data."
    in:
      reference_gtf_file: REFERENCE_GENOME_ANNOTATION
      seed_file: seed_file
      threads: THREADS
      #input_file_type: bamPairedEnd
      organism_name: organism_name
      metadata_file: metadata_file
      genome_seqs_dir: split_genome_seqs/output_dir # from seqkit process
      annotation_region_upstream: annotation_region_upstream
      annotation_region_downstream: annotation_region_downstream
      annotation_type: annotation_type
      bam_files: mapping_cageseq_data/aligned_bam_files # from star mapping process
      update_gtf_filename: update_gtf_filename
    out:
      - clustered_consensus_TSSs_bed
      - clustered_consensus_TSSs_txt
      - clustered_assigned_TSSs_txt
      - joined_clusters
      - all_cage_cluster_feature_uniq_gene_file
      - all_tss_feature_uniq_gene_file
      - all_tss_feature_file
      - output_gtf_file
      - stdout_log


#### outputs (upper-case) ####
outputs:
  # - id: SEQKIT_OUTPUT_DIR
  #   type: Directory
  #   label: "seqkit output directory"
  #   doc: "seqkit output directory"
  #   outputSource: split_genome_seqs/output_dir

  # - id: SPLIT_GENOME_SEQUENCE
  #   type: File[]
  #   label: "split genome sequence files"
  #   doc: "split genome sequence files"
  #   outputSource: split_genome_seqs/split_seqs


  # - id: TRIMMED_FASTQ1_FILES
  #   type: File[]
  #   label: "trimmed fastq1 files"
  #   doc: "trimmed fastq1 files"
  #   format: edam:format_1930 # fastq format
  #   outputSource: trimming_fastp_pe/trimmed_fastq1_files

  # - id: TRIMMED_FASTQ2_FILES
  #   type: File[]
  #   label: "trimmed fastq2 files"
  #   doc: "trimmed fastq2 files"
  #   format: edam:format_1930 # fastq format
  #   outputSource: trimming_fastp_pe/trimmed_fastq2_files

  - id: TRIMMING_REPORT_HTML_FILES
    type: File[]
    label: "trimming report html files"
    doc: "trimming report html files"
    format: edam:format_2331 # html report
    outputSource: trimming_fastp_pe/trimming_report_html_files

  - id: TRIMMING_REPORT_JSON_FILES
    type: File[]
    label: "trimming report json files"
    doc: "trimming report json files"
    format: edam:format_3464 # json
    outputSource: trimming_fastp_pe/trimming_report_json_files

  - id: FASTP_STDOUT_LOG_FILES
    type: File[]
    label: "fastp stdout log files"
    doc: "fastp stdout log files"
    format: edam:format_3671 # text
    outputSource: trimming_fastp_pe/stdout_log_files

  - id: FASTP_STDERR_LOG_FILES
    type: File[]
    label: "fastp stderr log files"
    doc: "fastp stderr log files"
    format: edam:format_3671 # text
    outputSource: trimming_fastp_pe/stderr_log_files


  # - id: STAR_INDEX_DIR
  #   type: Directory
  #   label: "STAR index directory"
  #   doc: "STAR index directory"
  #   outputSource: create_star_index/star_index_dir

  - id: STAR_INDEX_STDOUT_LOG
    type: File
    label: "STAR index stdout log"
    doc: "STAR index stdout log"
    format: edam:format_3671 # text
    outputSource: create_star_index/stdout_log

  - id: STAR_INDEX_STDERR_LOG
    type: File
    label: "STAR index stderr log"
    doc: "STAR index stderr log"
    format: edam:format_3671 # text
    outputSource: create_star_index/stderr_log

  - id: ALIGNED_BAM_FILES
    type: File[]
    label: "aligned BAM files"
    doc: "aligned BAM files"
    format: edam:format_2572 # BAM format
    outputSource: mapping_cageseq_data/aligned_bam_files

  - id: STAR_MAPPING_FINAL_LOG_FILES
    type: File[]
    label: "STAR mapping final log files"
    doc: "STAR mapping final log files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/final_log_files

  - id: STAR_MAPPING_MAIN_LOG_FILES
    type: File[]
    label: "STAR mapping main log files"
    doc: "STAR mapping main log files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/main_log_files


  - id: STAR_MAPPING_PROGRESS_LOG_FILES
    type: File[]
    label: "STAR mapping progress log files"
    doc: "STAR mapping progress log files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/progress_log_files


  - id: STAR_MAPPING_SJ_TAB_FILES
    type: File[]
    label: "STAR mapping splice junctions tab files"
    doc: "STAR mapping splice junctions tab files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/sj_tab_files

  - id: STAR_MAPPING_LOG_STDOUT_FILES
    type: File[]
    label: "STAR mapping stdout log files"
    doc: "STAR mapping stdout log files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/log_stdout_files

  - id: STAR_MAPPING_LOG_STDERR_FILES
    type: File[]
    label: "STAR mapping stderr log files"
    doc: "STAR mapping stderr log files"
    format: edam:format_3671 # text
    outputSource: mapping_cageseq_data/log_stderr_files


  - id: CLUSTERED_CONSENSUS_TSSs_BED
    type: File[]
    label: "clustered consensus TSSs bed file"
    doc: "clustered consensus TSSs bed file"
    format: edam:format_2572 # bed file
    outputSource: update_gtf/clustered_consensus_TSSs_bed

  - id: CLUSTERED_CONSENSUS_TSSs_TXT
    type: File[]
    label: "clustered consensus TSSs txt file"
    doc: "clustered consensus TSSs txt file"
    format: edam:format_3671 # text
    outputSource: update_gtf/clustered_consensus_TSSs_txt
    
  - id: CLUSTERED_ASSIGNED_TSSs_TXT
    type: File[]
    label: "clustered assigned TSSs txt file"
    doc: "clustered assigned TSSs txt file"
    format: edam:format_3671 # text
    outputSource: update_gtf/clustered_assigned_TSSs_txt
    
  - id: JOINED_CLUSTERS
    type: File
    label: "joined assigned clusters"
    doc: "joined assigned clusters file containing merged data from all input files"
    format: edam:format_3475 # tsv
    outputSource: update_gtf/joined_clusters
    
  - id: ALL_CAGE_CLUSTER_FEATURE_UNIQ_GENE_FILE
    type: File
    label: "all cage cluster feature uniq gene file"
    doc: "Contains unique cluster information for each gene, extracted from all_tss_feature.tsv."
    format: edam:format_3475 # tsv
    outputSource: update_gtf/all_cage_cluster_feature_uniq_gene_file

  - id: ALL_TSS_FEATURE_UNIQ_GENE_FILE
    type: File
    label: "all tss feature uniq gene file"
    doc: "Contains unique TSS information for each gene, extracted from all_tss_feature.tsv"
    format: edam:format_3475 # tsv
    outputSource: update_gtf/all_tss_feature_uniq_gene_file
    
  - id: ALL_TSS_FEATURE_FILE
    type: File
    label: "all tss feature uniq gene file"
    doc: "Contains TSS cluster information filtered based on the tag accumulation at the dominant TSS (tags.dominant_tss). Includes clusters from groups with the highest tag accumulation."
    format: edam:format_3475 # tsv
    outputSource: update_gtf/all_tss_feature_file

  - id: OUTPUT_GTF_FILE
    type: File
    label: "updated GFF format file"
    doc: "updated GFF format file with TSS (transcription start sites) information added"
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    outputSource: update_gtf/output_gtf_file

  - id: UPDATE_GTF_STDOUT_LOG
    type: File
    label: "update gtf stdout log"
    doc: "update gtf stdout log"
    outputSource: update_gtf/stdout_log


$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/