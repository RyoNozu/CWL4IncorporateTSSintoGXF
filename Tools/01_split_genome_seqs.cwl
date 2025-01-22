#!/usr/bin/env cwl-runner
# Generated from: seqkit ./Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz -i -O ./seqs_by_scaffold
class: CommandLineTool
cwlVersion: v1.2
label: "split genome seqs by scaffold"
doc: split genome seqs by scaffold. using seqkit version 2.8.2.

requirements:
  ShellCommandRequirement: {}


inputs:
  - id: genome_seqs
    type: File
    format: edam:format_3989
    label: "genome sequence gzipped file"
    doc: "genome sequence gzipped file"
    default:
      class: File
      format: edam:format_3989
      location: ../Data/Halichoeres_trimaculatus/Halichoeres_trimaculatus-hifiasm-3ddna-v1.1.edit.fna.gz
  - id: output_dir_name
    type: string
    label: "output directory name"
    doc: "output directory name"
    default: "seqs_by_scaffold"


arguments:
  - shellQuote: false
    valueFrom: |
      mkdir -p $(inputs.output_dir_name) && seqkit split $(inputs.genome_seqs.path) -i -O $(inputs.output_dir_name)


outputs:
  - id: all-for-debugging
    type:
      type: array
      items: [File, Directory]
    outputBinding:
      glob: "*"


hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
