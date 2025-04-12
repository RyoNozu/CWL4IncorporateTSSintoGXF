#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "make star index"
doc: "create star index for reference genome (STAR version 2.7.11b)"
requirements:
  ShellCommandRequirement: {}

inputs:
  - id: output_dir_name
    type: string
    label: "output directory name"
    doc: "output directory name"
    default: "star_genome_idx"
  - id: reference_genome
    type: File
    label: "reference genome fasta file"
    doc: "reference genome fasta file"
    format: edam:format_1929
    default:
      class: File
      format: edam:format_1929
      location: ../Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna
  - id: reference_genome_annotation
    type: File
    label: "reference genome annotation file"
    doc: "reference genome annotation file"
    format: edam:format_2306
    default:
      class: File
      format: edam:format_2306
      location: ../Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf
  - id: star_threads
    type: int
    label: "threads for star"
    doc: "threads for star"
    default: 16

stdout: prepare_star-index.log
stderr: prepare_star-index.log


arguments:
  - shellQuote: false
    valueFrom: |
      echo "Start to make STAR index for Htri reference genome gtf file" "`date '+%Y/%m/%d %H:%M:%S'`"
      STAR --version
      mkdir --verbose --parents $(inputs.output_dir_name)
      STAR --runMode genomeGenerate \
      --genomeDir $(inputs.output_dir_name) \
      --genomeFastaFiles $(inputs.reference_genome.path) \
      --sjdbGTFfile $(inputs.reference_genome_annotation.path) \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbOverhang 100 \
      --runThreadN $(inputs.star_threads)

outputs:
  - id: index_output_dir
    type: Directory
    label: "STAR index output directory"
    doc: "STAR index output directory (including text files, log files, etc.)"
    outputBinding:
      glob: $(inputs.output_dir_name)
  - id: stdout_log
    type: File
    label: "stdout log"
    doc: "stdout log"
    outputBinding:
      glob: prepare_star-index.log
  - id: stderr_log
    type: File
    label: "stderr log"
    doc: "stderr log"
    outputBinding:
      glob: prepare_star-index.log
      
      

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/star:2.7.11b--h5ca1c30_5

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/