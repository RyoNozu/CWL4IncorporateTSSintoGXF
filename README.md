# CWL4IncorporateTSSintoGXF

-   

## Requirements

- [cwltool](https://github.com/common-workflow-language/cwltool)    

```
# for example
conda create -n cwltool
conda activate cwltool
conda install -c conda-forge cwltool
```

- docker

## Simple usage
- clone this repository  
```
git clone https://github.com/RyoNozu/CWL4IncorporateTSSintoGXF.git
cd CWL4IncorporateTSSintoGXF
```
- run workflow  
```
cwltool --debug --outdir ./out/ ./workflow/cageseq_gtf_update_pe.cwl ./config/Commandlinetool_config/06_tssr_config.yaml
```

## Input files  

- CAGE-seq Read (fastq, paried/single-end)  
- reference genome (fasta)  
- gene annotation file (gff/gtf)  
- (BSgenome_data_package_seed_file (.txt))  
        > refere to forgeBSgenomeDataPkg function in [BSgenomeForge](https://bioconductor.org/packages/release/bioc/html/BSgenomeForge.html) packag  

## Output files  

- updated gxf file (.gff/gtf)  

***
