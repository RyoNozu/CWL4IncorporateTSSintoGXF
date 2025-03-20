#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "create STAR index"
doc: "create STAR index for mapping CAGE-Seq data (step 1: decompress reference genome fasta file, step 2: create STAR index)"

inputs:
  - id: create_star_index_threads
    type: int
    label: "threads (pigz and STAR command)"
    doc: "threads for decompressing reference genome fasta file and creating STAR index"
    default: 4

  - id: input_compressed_file
    type: File
    label: "compressed fasta file"
    doc: "input compressed reference genome fasta file (gzipped)"
    format: edam:format_3989
    default:
      class: File
      format: edam:format_3989
      path: ../Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz


  - id: star_index_dir_name
    type: string
    label: "output directory name"
    doc: "output directory name"
    default: "star_genome_idx"
  
  - id: reference_genome_annotation
    type: File
    label: "reference genome annotation file"
    doc: "reference genome annotation file"
    format: edam:format_2306
    default:
      class: File
      format: edam:format_2306
      path: ../Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf


steps:
  - id: decompress_reference_genome
    run: ../Tools/03_decompressed_pigz.cwl
    in:
      pigz_threads: create_star_index_threads
      input_compressed_file: input_compressed_file
    out: [output_decompressed_file]

  - id: create_star_index
    run: ../Tools/04_make_star_index.cwl
    in:
      output_dir_name: star_index_dir_name
      reference_genome: decompress_reference_genome/output_decompressed_file
      reference_genome_annotation: reference_genome_annotation
      star_threads: create_star_index_threads
    out: [index_output_dir, stdout_log, stderr_log]

outputs:
  - id: star_index_dir
    type: Directory
    label: "output directory for STAR index"
    doc: "output directory for STAR index"
    outputSource: create_star_index/index_output_dir

  - id: stdout_log
    type: File
    label: "stdout log"
    doc: "stdout log"
    format: edam:format_3671
    outputSource: create_star_index/stdout_log

  - id: stderr_log
    type: File
    label: "stderr log"
    doc: "stderr log"
    format: edam:format_3671
    outputSource: create_star_index/stderr_log

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/