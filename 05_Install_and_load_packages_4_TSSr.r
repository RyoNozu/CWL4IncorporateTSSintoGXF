# Install & load packages for TSSr

install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("Gviz")
BiocManager::install("rtracklayer")
BiocManager::install("DESeq2")
BiocManager::install("BSgenome")

install.packages("data.table")
install.packages("stringr")

devtools::install_github("Linlab-slu/TSSr", build_vignettes = TRUE)

library(TSSr)