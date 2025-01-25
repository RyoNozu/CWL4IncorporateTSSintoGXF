#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "trimming by fastp"
doc: "trimming by fastp (for pair-end CAGE-seq fastq data). using fastp version 0.23.2. Modified from https://github.com/nigyta/bact_genome/blob/master/cwl/tool/fastp/fastp.cwl"
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
      label: "fastq1"
      default:
        class: File
        path: ../Data/Halichoeres_trimaculatus/CAGE/All_data/Fastq/MK.F1/MK.F1_R1.fastq.gz
    - id: fastq2
      type: File
      label: "fastq2"
      default:
        class: File
        path: ../Data/Halichoeres_trimaculatus/CAGE/All_data/Fastq/MK.F1/MK.F1_R2.fastq.gz
    - id: compression_level
      type: int
      default: 9
    - id: threads
      type: int
      default: 16

outputs:
    - id: out_fastq1
      type: File
      label: "trimmed fastq1 file"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - id: out_fastq2
      type: File
      label: "trimmed fastq2 file"
      outputBinding:
        glob: $(inputs.fastq2.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.fq.gz
    - id: out_html
      type: File
      format: edam:format_2331
      label: "trimming report html file"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.html
    - id: out_json
      type: File
      format: edam:format_3464
      label: "trimming report json file"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_trim.report.json
    - id: stderr_log
      type: File
      format: edam:data_3671
      label: "stderr log file"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stderr_run_fastp.log
    - id: stdout_log
      type: File
      format: edam:data_3671
      label: "stdout log file"
      outputBinding:
        glob: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stdout_run_fastp.log

stderr: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stderr_run_fastp.log
stdout: $(inputs.fastq1.basename.replace(/\.gz$|\.bz2$/, '').replace(/\.fq$|\.fastq$/, ''))_stdout_run_fastp.log

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/fastp:0.23.2--hadf994f_5

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/