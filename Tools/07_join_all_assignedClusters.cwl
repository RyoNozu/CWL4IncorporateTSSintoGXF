#!/usr/bin/env cwl-runner
# Generated from: python3 09-01_join_all_assignedClusters.py *.assignedClusters.txt
class: CommandLineTool
cwlVersion: v1.2
label: "join all assignedClusters"
doc: "This custom python script combines all assigned Clusters.txt files processed by TSSr."

baseCommand: [python3]

inputs:
  - id: process_assigned_clusters_script
    type: File
    format: edam:format_3996 # python script
    label: "process assigned clusters custom python script"
    doc: "custom python script to process assigned clusters (TSSr output)"
    default:
      class: File
      format: edam:format_3996 # python script
      location: ../scripts/09-01_join_all_assignedClusters_modified.py
    inputBinding:
      position: 1

  - id: assigned_clusters_files
    type: File[]
    format: edam:format_3671 # text
    label: "assigned clusters files"
    doc: "assigned clusters files processed by TSSr"
    inputBinding:
      position: 2
      separate: true


outputs:
  - id: joined_clusters
    type: File
    label: "joined assigned clusters"
    doc: "joined assigned clusters file containing merged data from all input files"
    format: edam:format_3475 # tsv
    outputBinding:
      glob: "all-joined.assignedClusters.tsv"


hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pandas:2.2.1

$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/
