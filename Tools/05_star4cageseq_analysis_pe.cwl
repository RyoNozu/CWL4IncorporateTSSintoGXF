#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "run STAR for CAGE-Seq analysis (paired-end)"
doc: "Mapping the reads to the genome using STAR (STAR version 2.7.11b) for paired-end CAGE-seq data"

requirements:
  ShellCommandRequirement: {}

inputs:
  - id: star_index_dir
    type: Directory
    label: "STAR index directory"
    doc: "STAR index directory"

  - id: cage_seq_read_1
    type: File
    label: "CAGE-Seq read 1"
    doc: "CAGE-Seq read FASTQ file 1 (already trimmed)"
    format: edam:format_1930

  - id: cage_seq_read_2
    type: File
    label: "CAGE-Seq read 2"
    doc: "CAGE-Seq read FASTQ file 2 (already trimmed)"
    format: edam:format_1930

  - id: star_threads
    type: int
    label: "threads for star"
    doc: "threads for star"
    default: 16

stdout: "$(inputs.cage_seq_read_1.nameroot)_run_star.log"
stderr: "$(inputs.cage_seq_read_1.nameroot)_run_star.log"

 # Changed from zless command to zcat command for fastq.gz file
arguments:
  - shellQuote: false
    valueFrom: |
      echo "run STAR for Htri selected reference genome and gtf file" `date '+%Y/%m/%d %H:%M:%S'`
      STAR --version
      STAR --genomeDir $(inputs.star_index_dir.path)/ \
      --readFilesCommand zcat \
      --readFilesIn $(inputs.cage_seq_read_1.path) $(inputs.cage_seq_read_2.path) \
      --runThreadN $(inputs.star_threads) \
      --alignEndsType Extend5pOfRead1 \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix $(inputs.cage_seq_read_1.nameroot)_

outputs:
  - id: aligned_bam
    type: File
    label: "STAR aligned BAM file"
    doc: "STAR aligned BAM file"
    format: edam:format_2572  # BAM format
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_Aligned.sortedByCoord.out.bam"
  
  - id: final_log
    type: File
    label: "STAR final log file"
    doc: "STAR final log file"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_Log.final.out"
      
  - id: main_log
    type: File
    label: "STAR main log file"
    doc: "STAR main log file"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_Log.out"
      
  - id: progress_log
    type: File
    label: "STAR progress log file"
    doc: "STAR progress log file"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_Log.progress.out"
      
  - id: sj_tab
    type: File
    label: "STAR splice junctions file"
    doc: "STAR splice junctions file"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_SJ.out.tab"

  - id: log_stdout
    type: File
    label: "STAR stdout"
    doc: "STAR stdout"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_run_star.log"

  - id: log_stderr
    type: File
    label: "STAR stderr"
    doc: "STAR stderr"
    format: edam:format_3671
    outputBinding:
      glob: "$(inputs.cage_seq_read_1.nameroot)_run_star.log"


hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/star:2.7.11b--h5ca1c30_5

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/