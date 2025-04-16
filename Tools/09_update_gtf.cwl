#!/usr/bin/env cwl-runner
# Generated from: python3 10_update_gtf.py --gtf_file <input_gtf_file> --tss_file <input_tss_file> --output_file <output_gtf_file>
class: CommandLineTool
cwlVersion: v1.2
label: "update gtf files with TSS (Transcription Start Site) features"
doc: |
  This custom Python script to create a GTF/GFF file with additional TSS information.



baseCommand: [python3]
arguments:
  - $(inputs.update_gtf_script)
  - --gtf_file
  - $(inputs.gtf_file)
  - --tss_file
  - $(inputs.tss_file)
  - --output_file
  - $(inputs.update_gtf_filename)


inputs:
  - id: update_gtf_script
    type: File
    format: edam:format_3996 # python script
    label: "Custom Python script"
    doc: "Custom Python script to create a GTF/GFF file with additional TSS information"
    default:
      class: File
      format: edam:format_3996 # python script
      location: ../scripts/10_update_gtf.py

  - id: gtf_file
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "GTF/GFF file"
    doc: "GTF/GFF file"
    default:
      class: File
      format: edam:format_2305 # GFF format (including GTF, GFF3)
      location: ../Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf

  - id: tss_file
    type: File
    format: edam:format_3671 # text
    label: "TSS file"
    doc: "TSS file"
    default:
      class: File
      format: edam:format_3671 # text
      location: ../out/all_tss_feature_uniq.gene.tsv

  - id: update_gtf_filename
    type: string
    label: "File name in GFF format"
    doc: "File name in GFF format with TSS (transcription start sites) information added"
    default: "add_tss_feature.gtf"

stdout: update_gtf.log

outputs:
  - id: output_gtf_file
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "updated GFF format file"
    doc: "updated GFF format file with TSS (transcription start sites) information added"
    outputBinding:
      glob: $(inputs.update_gtf_filename)
  
  - id: stdout_log
    type: File
    label: "stdout log"
    doc: "stdout log"
    outputBinding:
      glob: update_gtf.log
      

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pandas:2.2.1


$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/