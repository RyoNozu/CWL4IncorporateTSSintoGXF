#!/usr/bin/env cwl-runner
# Generated from: python3 09-02_uniq_tss_feature_modified.py --input all-joined.assignedClusters.tsv
class: CommandLineTool
cwlVersion: v1.2
label: "uniq tss feature"
doc: |
  This custom python script extracts unique TSS features from all-joined.assignedClusters.tsv.
  Outputs are:
    - all_cage_cluster_feature_uniq.gene.tsv
    - all_tss_feature_uniq.gene.tsv
    - all_tss_feature.tsv

baseCommand: [python3]
arguments:
  - $(inputs.uniq_tss_feature_script)
  - --input
  - $(inputs.all_joined_assigned_clusters_file)

inputs:
  - id: uniq_tss_feature_script
    type: File
    format: edam:format_3996 # python script
    label: "uniq tss feature custom python script"
    doc: "custom python script to extract unique TSS features from all-joined.assignedClusters.tsv"
    default:
      class: File
      format: edam:format_3996 # python script
      location: ../scripts/09-02_uniq_tss_feature_modified_v2.py
    
  - id: all_joined_assigned_clusters_file
    type: File
    format: edam:format_3475 # tsv
    label: "all joined assigned clusters file"
    doc: "all joined assigned clusters file"

outputs:
  - id: all_cage_cluster_feature_uniq_gene_file
    type: File
    label: "all cage cluster feature uniq gene file"
    doc: "Contains unique cluster information for each gene, extracted from all_tss_feature.tsv."
    format: edam:format_3475 # tsv
    outputBinding:
      glob: "all_cage_cluster_feature_quniq.gene.tsv"

  - id: all_tss_feature_uniq_gene_file
    type: File
    label: "all tss feature uniq gene file"
    doc: "Contains unique TSS information for each gene, extracted from all_tss_feature.tsv"
    format: edam:format_3475 # tsv
    outputBinding:
      glob: "all_tss_feature_uniq.gene.tsv"

  - id: all_tss_feature_file
    type: File
    label: "all tss feature uniq gene file"
    doc: "Contains TSS cluster information filtered based on the tag accumulation at the dominant TSS (tags.dominant_tss). Includes clusters from groups with the highest tag accumulation."
    format: edam:format_3475 # tsv
    outputBinding:
      glob: "all_tss_feature.tsv"

hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pandas:2.2.1  
  
$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/
