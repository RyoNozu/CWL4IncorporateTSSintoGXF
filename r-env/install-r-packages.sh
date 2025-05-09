#!/bin/bash

# Install and load packages for TSSr execution
R -e "
# BiocManager setup (R 4.4.1)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager', version = '3.20')

# Bioconductor (v3.20) packages
BiocManager::install(c(\"Rsamtools\", \"GenomicRanges\", \"GenomicFeatures\", \"Gviz\", \"rtracklayer\", \"DESeq2\", \"BSgenome\", \"BSgenomeForge\", \"data.table\", \"stringr\", \"txdbmaker\"))

# Install from GitHub
# Install required packages with specific versions
install.packages('optparse', version = '1.7.5')
install.packages('devtools', version = '2.4.5')
devtools::install_github(\"Linlab-slu/TSSr\", ref = \"v0.99.1\", build_vignettes = TRUE, force = TRUE)

# Print installed package versions for verification
installed.packages()[c('devtools', 'BiocManager', 'Rsamtools', 'GenomicRanges', 
                       'GenomicFeatures', 'Gviz', 'rtracklayer', 'DESeq2', 
                       'BSgenome', 'BSgenomeForge', 'data.table', 'stringr', 'txdbmaker'), 'Version']
"