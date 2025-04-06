#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "trimming fastq files (paired-end)"
doc: "multiple fastq files trimming process using fastp version 0.23.4 and scatter feature requirement"

requirements:
  # https://www.commonwl.org/user_guide/topics/workflows.html#scattering-steps
  ScatterFeatureRequirement: {}

inputs:
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

    - id: compression_level
      type: int
      label: "compression level"
      doc: "compression level (default: 9)"
      default: 9
      
    - id: threads
      type: int
      label: "threads"
      doc: "threads (default: 16)"
      default: 16

steps:
    - id: trimming_fastp
      run: ../Tools/02_trimming_fastp_pe.cwl
      scatter: [fastq1, fastq2] # Parameters in Tool/02_trimming_fastp_pe.cwl should be listed here
      scatterMethod: dotproduct
      in:
        fastq1: fastq1_files
        fastq2: fastq2_files
        compression_level: compression_level
        threads: threads
      out: [out_fastq1, out_fastq2, out_html, out_json, stdout_log, stderr_log]

outputs:
    - id: trimmed_fastq1_files
      type: File[]
      format: edam:format_1930
      label: "trimmed fastq1 files"
      doc: "trimmed fastq1 files"
      outputSource: trimming_fastp/out_fastq1

    - id: trimmed_fastq2_files
      type: File[]
      format: edam:format_1930
      label: "trimmed fastq2 files"
      doc: "trimmed fastq2 files"
      outputSource: trimming_fastp/out_fastq2

    - id: trimming_report_html_files
      type: File[]
      format: edam:format_2331
      label: "trimming report html files"
      doc: "trimming report html files"
      outputSource: trimming_fastp/out_html

    - id: trimming_report_json_files
      type: File[]
      format: edam:format_3464
      label: "trimming report json files"
      doc: "trimming report json files"
      outputSource: trimming_fastp/out_json

    - id: stdout_log_files
      type: File[]
      format: edam:data_3671
      label: "stdout log files"
      doc: "stdout log files"
      outputSource: trimming_fastp/stdout_log 

    - id: stderr_log_files
      type: File[]
      format: edam:data_3671
      label: "stderr log files"
      doc: "stderr log files"
      outputSource: trimming_fastp/stderr_log

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/