#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "decompress reference genome fasta file"
doc: "decompress reference genome fasta file using pigz"

requirements:
  ShellCommandRequirement: {}

baseCommand: [pigz]
arguments:
  - --decompress
  - --keep
  - --stdout
  - --processes
  - $(inputs.pigz_threads)
  - $(inputs.input_compressed_file)


inputs:
  - id: pigz_threads
    type: int
    label: "pigz threads"
    doc: "pigz threads"
    default: 16

  - id: input_compressed_file
    type: File
    label: "input compressed reference genome fasta file"
    doc: "input compressed reference genome fasta file"
    format: edam:format_3989
    default:
      class: File
      format: edam:format_3989
      path: ../Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz

stdout: $(inputs.input_compressed_file.nameroot)


outputs:
  - id: output_decompressed_file
    type: File
    label: "output decompressed reference genome fasta file"
    doc: "output decompressed reference genome fasta file"
    format: edam:format_1929
    outputBinding:
      glob: $(inputs.input_compressed_file.nameroot)

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pigz:2.8

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/