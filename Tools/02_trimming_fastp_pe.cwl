#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "trimming by fastp"
doc: "trimming by fastp (for pair-end CAGE-seq fastq data). using fastp version 0.23.4. Modified from https://github.com/nigyta/bact_genome/blob/master/cwl/tool/fastp/fastp.cwl"
requirements:
  InlineJavascriptRequirement: {}

baseCommand: [fastp]
arguments:
    - -i
    - $(inputs.fastq1.path)
    - -I
    - $(inputs.fastq2.path)
    - -o
    - $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - -O
    - $(inputs.fastq2.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - -h
    - $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.html
    - -j
    - $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.json
    - -Q
    - -z
    - $(inputs.compression_level)
    - -w
    - $(inputs.threads)

inputs:
    - id: fastq1
      type: File
      format: edam:format_1930
      label: "fastq1"
      doc: "fastq1 file (e.g. MK.F1_R1.fastq.gz)"

    - id: fastq2
      type: File
      format: edam:format_1930
      label: "fastq2"
      doc: "fastq2 file (e.g. MK.F1_R2.fastq.gz)"

    - id: compression_level
      label: "compression level"
      doc: "compression level (default: 9)"
      type: int
      default: 9
    - id: threads
      label: "threads"
      doc: "threads (default: 16)"
      type: int
      default: 16

outputs:
    - id: out_fastq1
      type: File
      format: edam:format_1930
      label: "trimmed fastq1 file"
      doc: "trimmed paired-end fastq file (e.g. MK.F1_R1_trim.fq.gz)"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - id: out_fastq2
      type: File
      format: edam:format_1930
      label: "trimmed fastq2 file"
      doc: "trimmed paired-end fastq file (e.g. MK.F1_R2_trim.fq.gz)"
      outputBinding:
        glob: $(inputs.fastq2.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - id: out_html
      type: File
      format: edam:format_2331
      label: "trimming report html file"
      doc: "trimming report html file (e.g. MK.F1_R1_trim.report.html)"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.html
    - id: out_json
      type: File
      format: edam:format_3464
      label: "trimming report json file"
      doc: "trimming report json file (e.g. MK.F1_R1_trim.report.json)"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.json
    - id: stdout_log
      type: File
      format: edam:data_3671
      label: "stdout log file"
      doc: "stdout log file (e.g. MK.F1_R1_stdout_run_fastp.log)"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stdout_run_fastp.log
    - id: stderr_log
      type: File
      format: edam:data_3671
      label: "stderr log file"
      doc: "stderr log file (e.g. MK.F1_R1_stderr_run_fastp.log)"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stderr_run_fastp.log

stdout: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stdout_run_fastp.log
stderr: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stderr_run_fastp.log

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/fastp:0.23.4--h125f33a_4

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/