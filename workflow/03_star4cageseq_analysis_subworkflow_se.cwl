#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "step4: STAR for CAGE-Seq reads analysis"
doc: "STAR for CAGE-Seq reads analysis"

requirements:
  # https://www.commonwl.org/user_guide/topics/workflows.html#scattering-steps
  ScatterFeatureRequirement: {}

inputs:
  - id: star_index_dir
    type: Directory
    label: "STAR index directory"
    doc: "STAR index directory"
    default:
      class: Directory
      location: ../out/star_genome_idx

  - id: cage_seq_reads
    type: File[]
    format: edam:format_1930
    label: "CAGE-Seq reads"
    doc: "CAGE-Seq reads"

  - id: star_threads
    type: int
    label: "threads for star"
    doc: "threads for star"
    default: 16

steps:
  - id: star_analysis_se
    run: ../Tools/05_star4cageseq_analysis_se.cwl
    scatter: [cage_seq_read] # Parameters in Tool/05_star4cageseq_analysis_se.cwl should be listed here
    scatterMethod: dotproduct
    in:
      star_index_dir: star_index_dir
      cage_seq_read: cage_seq_reads
      star_threads: star_threads
    out: [aligned_bam, final_log, main_log, progress_log, sj_tab, log_stdout, log_stderr]

outputs:
  - id: aligned_bam_files
    type: File[]
    label: "STAR aligned BAM files"
    doc: "STAR aligned BAM files"
    format: edam:format_2572
    outputSource: star_analysis_se/aligned_bam

  - id: final_log_files
    type: File[]
    label: "STAR final log files"
    doc: "STAR final log files"
    format: edam:format_3671
    outputSource: star_analysis_se/final_log

  - id: main_log_files
    type: File[]
    label: "STAR main log files"
    doc: "STAR main log files"
    format: edam:format_3671
    outputSource: star_analysis_se/main_log

  - id: progress_log_files
    type: File[]
    label: "STAR progress log files"
    doc: "STAR progress log files"
    format: edam:format_3671
    outputSource: star_analysis_se/progress_log

  - id: sj_tab_files
    type: File[]
    label: "STAR splice junctions tab files"
    doc: "STAR splice junctions tab files"
    format: edam:format_3671
    outputSource: star_analysis_se/sj_tab

  - id: log_stdout_files
    type: File[]
    label: "STAR stdout files"
    doc: "STAR stdout files"
    format: edam:format_3671
    outputSource: star_analysis_se/log_stdout

  - id: log_stderr_files
    type: File[]
    label: "STAR stderr files"
    doc: "STAR stderr files"
    format: edam:format_3671
    outputSource: star_analysis_se/log_stderr

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/