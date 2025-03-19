#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: "mapping CAGEseq data and update gtf file"
doc: "mapping CAGEseq data and update gtf file"

requirements:
  # https://www.commonwl.org/user_guide/topics/workflows.html#scattering-steps
  ScatterFeatureRequirement: {}


inputs:


steps:
  - id: processing_fastq
    run: ./01_trimming_fastq_subworkflow_pe.cwl