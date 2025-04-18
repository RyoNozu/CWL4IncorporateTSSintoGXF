#!/usr/bin/env cwl-runner
# Generated from: Rscript ./scripts/05-08_combined_exec_TSSr_modified.R --refSource ./Data/Halichoeres_trimaculatus/braker_correctID_v3.gtf --seedFile ./Data/Halichoeres_trimaculatus/BSgenome.Htrimaculatus.Htriv1.1-seed_2.txt --threads 16 --inputFilesType bamPairedEnd --organismName Halichoeres trimaculatus --metadata ./Data/Halichoeres_trimaculatus/sample_metadata.csv --bam_files ./out/MK.F1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.F4_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M1.Mix_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M2_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M3_R1_trim.fq_Aligned.sortedByCoord.out.bam,./out/MK.M4_R1_trim.fq_Aligned.sortedByCoord.out.bam
class: CommandLineTool
cwlVersion: v1.2
label: "combined exec TSSr"
doc: |
  TSSr (https://github.com/Linlab-slu/TSSr/tree/v0.99.1) execution Rscript for detecting TSSs (transcription start sites) from bam files processed by STAR mapping process.
  In this version, Dr. Ryo Nozu modified the TSSr repository (https://github.com/RyoNozu/TSSr).

baseCommand: [Rscript]
arguments:
  - $(inputs.Rscript)
  - --refSource
  - $(inputs.referenceFile)
  - --seedFile
  - $(inputs.seedFile)
  - --threads
  - $(inputs.threads)
  - --inputFilesType
  - $(inputs.inputFilesType)
  - --organismName
  - $(inputs.organismName)
  - --metadata
  - $(inputs.metadata)
  - --genome_seqs_dir
  - $(inputs.genome_seqs_dir)
  - --upstream
  - $(inputs.annotation_region_upstream)
  - --downstream
  - $(inputs.annotation_region_downstream)
  - --annotationType
  - $(inputs.annotation_type)

inputs:
  - id: Rscript
    type: File
    format: edam:format_3999 # R script
    label: "custom R script"
    doc: "custom R script for executing TSSr. TSSr repository is modified by Dr. Ryo Nozu (https://github.com/RyoNozu/TSSr)"
    default: 
      class: File
      format: edam:format_3999 # R script
      location: ../scripts/05-08_combined_exec_TSSr_modified.R

  - id: referenceFile
    type: File
    format: edam:format_2305 # GFF format (including GTF, GFF3)
    label: "reference file (GTF or GFF format)"
    doc: "reference file (GTF or GFF format)"

  - id: seedFile
    type: File
    format: edam:format_3671 # text
    label: "seed file"
    doc: "seed file for BSgenome package"

  - id: threads
    type: int
    label: "number of threads"
    doc: "number of threads"
    default: 16

  - id: inputFilesType
    type: string
    label: "input files type"
    doc: "input BAM files type (bamPairedEnd or bam). For example, if the processed CAGE-seq files are paired-end files, select “bamPairedEnd”."
    default: "bamPairedEnd"

  - id: organismName
    type: string
    label: "organism name"
    doc: "organism name"
    default: "Halichoeres trimaculatus"

  - id: metadata
    type: File
    format: edam:format_3752 # csv
    label: "metadata file"
    doc: "metadata file for grouping bam files"

  - id: genome_seqs_dir
    type: Directory
    label: "genome sequences directory"
    doc: "genome sequences directory processed by seqkit process"


  - id: annotation_region_upstream
    type: int
    label: "annotation region upstream"
    doc: "The maximum allowable distance for assigning a dominant TSS (transcription start site) to its downstream gene, when the TSS is located upstream of 5'-end of the CDS or the the transcript. Default: 1000 bp"
    default: 1000

  - id: annotation_region_downstream
    type: int
    label: "annotation region downstream"
    doc: "The maximum allowable distance for assigning a TSS when the dominant TSS is located downstream of 5'-end of the CDS or the the transcript. Default: 0 (i.e., TSSs located downstream of the 5'-end of feature (CDS or transcript) are not assigned to the gene)"
    default: 0

  - id: annotation_type
    type: string
    label: "annotation type"
    doc: "annotation type (\"genes\" or \"transcripts\")"
    default: "genes"

  - id: bam_files
    type: File[]
    format: edam:format_2572 # bam file
    label: "bam files"
    doc: "bam files processed by STAR mapping process"
    inputBinding:
      prefix: --bam_files
      itemSeparator: ","
      separate: true

outputs:
  - id: clustered_consensus_TSSs_bed
    type: File[]
    format: edam:format_2572 # bed file
    label: "clustered consensus TSSs bed file"
    doc: "clustered consensus TSSs bed file"
    outputBinding:
      glob: "*.consensusClusters.bed"

  - id: clustered_consensus_TSSs_txt
    type: File[]
    format: edam:format_3671 # text
    label: "clustered consensus TSSs txt file"
    doc: "clustered consensus TSSs txt file"
    outputBinding:
      glob: "*.consensusClusters.txt"

  - id: clustered_assigned_TSSs_txt
    type: File[]
    format: edam:format_3671 # text
    label: "clustered assigned TSSs txt file"
    doc: "clustered assigned TSSs txt file"
    outputBinding:
      glob: "*.assignedClusters.txt"

hints:
  - class: DockerRequirement
    # original Docker image: https://hub.docker.com/r/sorayone56/tssr-r-env
    dockerPull: sorayone56/tssr-r-env:1.0.2


$namespaces:
  s: https://schema.org/
  edam: https://edamontology.org/