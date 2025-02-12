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

# generate BSgenome Package for target species
# This is example case, H. trimaculatus
# seed file: BSgenome.Htrimaculatus.Htriv1.1-seed.txt
BSgenomeForge::forgeBSgenomeDataPkg(".+-seed.txt", replace = TRUE)
devtools::build("BSgenome.Htrimaculatus.inhouse.Htriv1")
devtools::check_built("BSgenome.Htrimaculatus.inhouse.Htriv1_1.0.0.tar.gz")
devtools::install_local("BSgenome.Htrimaculatus.inhouse.Htriv1_1.0.0.tar.gz")


# Import required files & TSS calling
inputFiles <- list.files(path = "./", pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
inputFilesType <- "bamPairedEnd"  # set “inputFilesType” as “bamPairedEnd” for paired-end BAM files, and as "TSStable" if the input file is a TSS table 

args <- commandArgs(trailingOnly = TRUE)
refSource <- args[1] # コマンドライン引数からrefSourceを取得
prefix1 <- args[2]   # コマンドライン引数から1つ目のprefixを取得
prefix2 <- args[3]   # コマンドライン引数から2つ目のprefixを取得
prefix3 <- args[4]   # コマンドライン引数から3つ目のprefixを取得
prefix4 <- args[5]   # コマンドライン引数から4つ目のprefixを取得

# 各グループのファイル数を指定
group_sizes <- c(3, 3, 1, 1)

# グループごとにファイルを分ける
group1_files <- inputFiles[1:group_sizes[1]]
group2_files <- inputFiles[(group_sizes[1] + 1):(group_sizes[1] + group_sizes[2])]
group3_files <- inputFiles[(group_sizes[1] + group_sizes[2] + 1):(group_sizes[1] + group_sizes[2] + group_sizes[3])]
group4_files <- inputFiles[(group_sizes[1] + group_sizes[2] + group_sizes[3] + 1):length(inputFiles)]

# サンプルラベルを生成
sampleLabels_group1 <- paste0(prefix1, sprintf("%02d", 1:length(group1_files)))
sampleLabels_group2 <- paste0(prefix2, sprintf("%02d", 1:length(group2_files)))
sampleLabels_group3 <- paste0(prefix3, sprintf("%02d", 1:length(group3_files)))
sampleLabels_group4 <- paste0(prefix4, sprintf("%02d", 1:length(group4_files)))
sampleLabels <- c(sampleLabels_group1, sampleLabels_group2, sampleLabels_group3, sampleLabels_group4)

myTSSr <- new("TSSr", genomeName = "BSgenome.Htrimaculatus.inhouse.Htriv1"
                 ,inputFiles = inputFiles
                 ,inputFilesType = inputFilesType
                 ,sampleLabels = sampleLabels
                 ,sampleLabelsMerged = c("group1", "group2", "group3", "group4")
                 ,mergeIndex = c(rep(1, length(group1_files)), rep(2, length(group2_files)), rep(3, length(group3_files)), rep(4, length(group4_files)))
                 ,refSource = refSource
                 ,organismName = "Halichoeres trimaculatus")


### TSS calling
## options
# sequencingQualityThreshold: minimum sequencing quality score for a base to be considered (default:10, specified in integer)
# mappingQualityThreshold: minimum mapping quality score for a read to be considered (default:20, specified in integer)
# softClippingAllowed: (default: false), Note: If you ran STAR with the default alignment parameters (without using --alignEndsType Extend5pOfRead1), 
#                                              please ensure that you set softclippingAllowed = TRUE when running TSSr. 
##

getTSS(myTSSr, sequencingQualityThreshold=20, mappingQualityThreshold=20)


# TSS clustering, consensus TSS identification and annotation

### Merging samples (biological replicates)
# Users can merge multiple samples (e.g., biological replicates) into previously defined groups with 
# mergeSamples function. The "mergeIndex" argument directs which samples will be merged and how the final 
# dataset will be ordered accordingly.
### Filter raw TSS counts
## options
# method: poisson: 
# method: TPM: 
#          tpmLow: Default value is 0.1
##

mergeSamples(myTSSr)
filterTSS(myTSSr, method = "poisson")
filterTSS(myTSSr, method = "TPM", tpmLow = 0.05) 

### Clustering TSSs to infer core promoters with “clusterTSS” function
# The “clusterTSS” function was designed to group neighboring TSSs into distinct TSS clusters (TCs), representing putative core promoters
## options
# method: peakclu:
# peakDistance: 
# extensionDistance:
# localThreshold:
# clusterThreshold:
# useMultiCore:
# numCores:
##

clusterTSS(myTSSr, method = "peakclu", peakDistance = 100, extensionDistance = 20, 
           localThreshold = 0.02, clusterThreshold = 1, 
           useMultiCore = TRUE, numCores = 10)


### Aggregating consensus TSS clusters with “consensusCluster” function.
# TSSr infers a set of consensus core promoters using the “consensusCluster” function to assign the same ID for TCs belong to the same core promoter
## options
# dis: (default: 50) Minimum distance between two peaks to be aggregated together into the same consensus cluster.
# useMultiCore: Logical indicating whether multiple cores are used or not. Default is FALSE.
# numCores: Number of cores. Default is NULL
##

consensusCluster(myTSSr, dis = 100, useMultiCore = TRUE, numCores = 10)

### export the consensus TCs to txt file
#
## options
# data: "tagClusters", "consensusClusters", "assigned", "unassigned". Default is "assigned"
##

exportClustersTable(myTSSr, data = "consensusClusters")
exportClustersToBed(myTSSr, data = "consensusClusters")	# consensusCluster()後の出力ではエラーだが、一通り実行後に再度行うと出力できた(20241026)

### Assigning TCs to downstream genes
# Assigning TCs to downstream genes as their core promoters is required for annotation of the 5’ boundaries of genomic features. 
# This process is also a prerequisite for further interrogations of regulated transcription initiation at the gene level. 
# TSSr offers the “annotateCluster” function to assign TCs to their downstream genes. The assignment of a TC to a gene is 
# based on the distance between the position of the dominant TSS of a TC and the annotated 5’ends of coding sequences (CDS) or transcripts.
### ### ver: analyze_poisson-TPM005-data_peakDis100_extDis20_localTh002_consDis100_assi5k.RData
## options
# filterCluster: logical, whether to filter out TCs with a low expression level (default: TRUE)
# filterClusterThreshold: numeric, the threshold of the expression level for filtering TCs (default: 0.02)
# annotationType: character, the type of genomic features to be annotated (default: "genes")
# upstream: numeric, the distance upstream of the TSS to be considered for annotation (default: 1000)
# upstreamOverlap: Upstream distance to the start position of annotation feature if overlapped with the upstream neighboring feature. Default value = 500. 
# downstream: Downstream distance to the start position of annotation feature. Default value = 0. 
#                   Note: if annotationType == "transctipt" or the gene annotations start from transcription start sites (TSSs), the recommended value = 500.
##

annotateCluster(myTSSr, clusters = "consensusClusters", filterCluster = TRUE, 
                filterClusterThreshold = 0.02, annotationType = "genes", 
                upstream = 5000, upstreamOverlap = 500, downstream = 0)

### Export the annotated TCs in each group to txt files
# exported files prefix is each "sampleLabels" (refer to 07_ R script)
# e.g. "sampleLabels".assignedClusters.txt

exportClustersTable(myTSSr, data = "assigned")