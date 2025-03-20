#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "mapping CAGEseq data and update gtf file"
doc: "mapping CAGEseq data and update gtf file"

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
    format: edam:format_3989
    label: "genome sequence file (gzipped fasta format) for mapping"
    doc: "genome sequence file (gzipped fasta format) for mapping"
    default:
      class: File
      format: edam:format_3989
      path: ../Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz


  #### seqkit inputs (lower-case) ####
  - id: seqkit_output_dir_name
    type: string
    label: "seqkit output directory name"
    doc: "seqkit output directory name"
    default: "seqs_by_scaffold"

  #### fastp inputs (paired-end) (lower-case) ####
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

  - id: fastp_compression_level
    type: int
    label: "fastp compression level"
    doc: "fastp compression level"
    default: 9

  #### STAR index inputs (lower-case) ####
  - id: star_index_dir_name
    type: string
    label: "STAR index directory name"
    doc: "STAR index directory name"
    default: "star_genome_idx"

  - id: reference_genome_annotation
    type: File
    format: edam:format_2306
    label: "reference genome annotation file"
    doc: "reference genome annotation file"
    default:
      class: File
      format: edam:format_2306
      path: ../Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf


steps:
  #### 1. split genome sequences by scaffold ####
  - id: split_genome_seqs
    run: ../Tools/01_split_genome_seqs.cwl
    in:
      genome_seqs: GENOME_FILE
      output_dir_name: seqkit_output_dir_name
    out: [output_dir, split_seqs]

  #### 2. trimming fastq files ####
  - id: trimming_fastp_pe
    run: ./01_trimming_fastq_subworkflow_pe.cwl
    in:
      fastq1_files: fastq1_files
      fastq2_files: fastq2_files
      compression_level: fastp_compression_level
      threads: THREADS
    out: [trimmed_fastq1_files, trimmed_fastq2_files, trimming_report_html_files, trimming_report_json_files, stdout_log_files, stderr_log_files]
  
  #### 3. create STAR index ####
  - id: create_star_index
    run: ./02_star_index_subworkflow.cwl
    in:
      create_star_index_threads: THREADS
      input_compressed_file: GENOME_FILE
      star_index_dir_name: star_index_dir_name
      reference_genome_annotation: reference_genome_annotation
    out: [star_index_dir, stdout_log, stderr_log]

  #### 4. mapping CAGEseq data####
  - id: mapping_cageseq_data
    run: ./03_star4cageseq_analysis_subworkflow_pe.cwl
    in:
      star_index_dir: create_star_index/star_index_dir
      cage_seq_reads_1: trimming_fastp_pe/trimmed_fastq1_files
      cage_seq_reads_2: trimming_fastp_pe/trimmed_fastq2_files
      star_threads: THREADS
    out: [aligned_bam_files, final_log_files, main_log_files, progress_log_files, sj_tab_files, log_stdout_files, log_stderr_files]


outputs:
  - id: seqkit_output_dir
    type: Directory
    label: "seqkit output directory"
    doc: "seqkit output directory"
    outputSource: split_genome_seqs/output_dir

  - id: split_genome_sequence
    type: File[]
    label: "split genome sequence files"
    doc: "split genome sequence files"
    outputSource: split_genome_seqs/split_seqs

  - id: trimmed_fastq1_files
    type: File[]
    label: "trimmed fastq1 files"
    doc: "trimmed fastq1 files"
    outputSource: trimming_fastp_pe/trimmed_fastq1_files


  - id: trimmed_fastq2_files
    type: File[]
    label: "trimmed fastq2 files"
    doc: "trimmed fastq2 files"
    outputSource: trimming_fastp_pe/trimmed_fastq2_files

  - id: trimming_report_html_files
    type: File[]
    label: "trimming report html files"
    doc: "trimming report html files"
    outputSource: trimming_fastp_pe/trimming_report_html_files

  - id: trimming_report_json_files
    type: File[]
    label: "trimming report json files"
    doc: "trimming report json files"
    outputSource: trimming_fastp_pe/trimming_report_json_files

  - id: fastp_stdout_log_files
    type: File[]
    label: "fastp stdout log files"
    doc: "fastp stdout log files"
    outputSource: trimming_fastp_pe/stdout_log_files

  - id: fastp_stderr_log_files
    type: File[]
    label: "fastp stderr log files"
    doc: "fastp stderr log files"
    outputSource: trimming_fastp_pe/stderr_log_files
    
  - id: star_index_dir
    type: Directory
    label: "STAR index directory"
    doc: "STAR index directory"
    outputSource: create_star_index/star_index_dir

  - id: star_index_stdout_log
    type: File
    label: "STAR index stdout log"
    doc: "STAR index stdout log"
    outputSource: create_star_index/stdout_log

  - id: star_index_stderr_log
    type: File
    label: "STAR index stderr log"
    doc: "STAR index stderr log"
    outputSource: create_star_index/stderr_log

  - id: aligned_bam_files
    type: File[]
    label: "aligned BAM files"
    doc: "aligned BAM files"
    outputSource: mapping_cageseq_data/aligned_bam_files

  - id: star_mapping_final_log_files
    type: File[]
    label: "STAR mapping final log files"
    doc: "STAR mapping final log files"
    outputSource: mapping_cageseq_data/final_log_files

  - id: star_mapping_main_log_files
    type: File[]
    label: "STAR mapping main log files"
    doc: "STAR mapping main log files"
    outputSource: mapping_cageseq_data/main_log_files


  - id: star_mapping_progress_log_files
    type: File[]
    label: "STAR mapping progress log files"
    doc: "STAR mapping progress log files"
    outputSource: mapping_cageseq_data/progress_log_files


  - id: star_mapping_sj_tab_files
    type: File[]
    label: "STAR mapping splice junctions tab files"
    doc: "STAR mapping splice junctions tab files"
    outputSource: mapping_cageseq_data/sj_tab_files

  - id: star_mapping_log_stdout_files
    type: File[]
    label: "STAR mapping stdout log files"
    doc: "STAR mapping stdout log files"
    outputSource: mapping_cageseq_data/log_stdout_files

  - id: star_mapping_log_stderr_files
    type: File[]
    label: "STAR mapping stderr log files"
    doc: "STAR mapping stderr log files"
    outputSource: mapping_cageseq_data/log_stderr_files


$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/