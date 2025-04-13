#!/usr/bin/env Rscript

# Install & load packages for TSSr

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager") # version 1.30.25
# BiocManager::install("Rsamtools") # version 2.20.0
# BiocManager::install("GenomicRanges") # version 1.56.2
# BiocManager::install("GenomicFeatures") # version 1.56.0
# BiocManager::install("Gviz") # version 1.48.0
# BiocManager::install("rtracklayer") # version 1.64.0
# BiocManager::install("DESeq2") # version 1.44.0
# BiocManager::install("BSgenome") # version 1.72.0
# BiocManager::install("BSgenomeForge")
# BiocManager::install("txdbmaker")
# install.packages("data.table") # version 1.16.4
# install.packages("stringr") # version 1.5.1
# install.packages('optparse', version = '1.7.5')
# install.packages('devtools', version = '2.4.5')

# devtools::install_github("Linlab-slu/TSSr", ref = "v0.99.1", build_vignettes = TRUE,force = TRUE) # version 0.99.6

library(TSSr)
library(optparse)

# Option list
option_list <- list(
  make_option(
    c("-r", "--refSource"),
    type = "character",
    default = NULL,
    help="Reference source file (GTF or GFF file)", 
    metavar="FILE"
  ),
  make_option(
    c("-s", "--seedFile"),
    type = "character",
    default = NULL,
    help = "Seed file for BSgenome",
    metavar = "FILE"
  ),
  make_option(
    c("-d", "--bam_dir"),
    type = "character",
    default = "./",
    help = "Directory containing BAM files (used if --bam_files is not specified)",
    metavar = "DIRECTORY"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 8,
    help = "Number of threads",
    metavar = "INTEGER"
  ),
  make_option(
    c("-y", "--inputFilesType"),
    type = "character",
    default = "bamPairedEnd",
    help = "Input files type, either 'bamPairedEnd' or 'bam'",
    metavar = "STRING"
  ),
  make_option(
    c("-o", "--organismName"),
    type = "character",
    default = NULL,
    help = "Full scientific name of the organism (e.g., 'Halichoeres trimaculatus')",
    metavar = "NAME"
  ),
  make_option(
    c("-b", "--bam_files"),
    type = "character",
    default = NULL,
    help = "Comma-separated list of BAM files. such as 'sample1.bam,sample2.bam,sample3.bam'",
    metavar = "FILES"
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Path to CSV metadata file defining sample groups (required)",
    metavar = "FILE"
  )
)

# Parse options, including positional arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = TRUE)

# Split parsed options and remaining arguments
parsed_options <- opt$options
remaining_args <- opt$args

# Process for BAM files
if (!is.null(parsed_options$bam_files) && parsed_options$bam_files != "") {
  # カンマ区切りのBAMファイルリストを分割
  inputFiles <- unlist(strsplit(parsed_options$bam_files, ","))
  
  # 存在確認とフルパス変換
  inputFiles <- sapply(inputFiles, function(f) {
    if (!file.exists(f) && !startsWith(f, "/")) {
      f <- file.path(getwd(), f)
    }
    if (!file.exists(f)) {
      warning(paste("File not found:", f))
    }
    return(f)
  })
} else {
  # Default: Search for files from the directory
  inputFiles <- list.files(path = parsed_options$bam_dir, pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
}

# Check existence
if (length(inputFiles) == 0) {
  stop("BAM files not found. Please specify BAM files after --bam_files or check --bam_dir option.")
}

# 入力ファイル一覧を表示
cat("Input BAM files found:\n")
for (i in 1:length(inputFiles)) {
  cat(sprintf("%d: %s\n", i, basename(inputFiles[i])))
}

# シードファイルからパッケージ名を直接抽出する
extract_package_name_from_seed <- function(seed_file) {
  lines <- readLines(seed_file, n = 1)  # 最初の行だけ読む
  if (length(lines) > 0 && grepl("^Package:", lines[1])) {
    # "Package: XXX" の形式から XXX を抽出
    return(gsub("^Package:\\s*", "", lines[1]))
  }
  return(NULL)
}

# シードファイルからパッケージ名を取得
seed_file <- parsed_options$seedFile
package_name <- extract_package_name_from_seed(seed_file)
# パッケージ名から生物種名（organism name）を抽出
organism_name <- sub("^BSgenome\\.([^\\.]+).*$", "\\1", package_name)

if (is.null(package_name)) {
  # パッケージ名が見つからない場合はファイル名から作成
  seed_basename <- tools::file_path_sans_ext(basename(seed_file))
  package_name <- paste0("BSgenome.", seed_basename)
  organism_name <- seed_basename
}

# BSgenomeパッケージを生成
BSgenomeForge::forgeBSgenomeDataPkg(seed_file, replace = TRUE)

# パッケージをビルドしてインストール
if (dir.exists(package_name)) {
  # パッケージディレクトリが存在する場合はビルド
  devtools::build(package_name)
  # tarballファイル名を取得
  tarball_files <- list.files(pattern = ".*\\.tar\\.gz$")
  tarball_file <- grep(package_name, tarball_files, value = TRUE)[1]
  if (is.na(tarball_file)) {
    tarball_file <- paste0(package_name, "_1.0.0.tar.gz")
  }
  devtools::install_local(tarball_file)
  
  # TSSrオブジェクト生成に使用するゲノム名
  genome_name <- package_name
} else {
  stop(paste("BSgenome package directory not found:", package_name))
}

# Import required files & TSS calling
inputFilesType <- parsed_options$inputFilesType
refSource <- parsed_options$refSource

# メタデータファイルが指定されていなければエラー
if (is.null(parsed_options$metadata)) {
  stop("Error: Metadata file (--metadata) is required. Please provide a CSV file with sample_id and prefix columns.")
}

# メタデータファイルが存在するか確認
if (!file.exists(parsed_options$metadata)) {
  stop(paste("Error: Metadata file not found:", parsed_options$metadata))
}

# メタデータファイルを読み込む
cat("Using metadata file for grouping:", parsed_options$metadata, "\n")
metadata <- read.csv(parsed_options$metadata, stringsAsFactors = FALSE)

# メタデータファイルの形式を検証
if (!all(c("sample_id", "prefix") %in% colnames(metadata))) {
  stop("Error: Metadata file must contain 'sample_id' and 'prefix' columns.")
}

# プレフィックスの一覧を取得（重複を削除）
unique_prefixes <- unique(metadata$prefix)
group_count <- length(unique_prefixes)
prefixes <- unique_prefixes

cat(sprintf("Detected %d groups from metadata:\n", group_count))
for (i in 1:group_count) {
  cat(sprintf("  Group %d: %s\n", i, prefixes[i]))
}

# BAMファイルをグループに分類（空のグループも保持）
group_files <- vector("list", group_count)
names(group_files) <- prefixes

# BAMファイルの使用状況を追跡
used_bams <- character(0)

for (i in 1:group_count) {
  # 現在のグループに属するサンプルIDを取得
  current_prefix <- prefixes[i]
  samples_in_group <- metadata$sample_id[metadata$prefix == current_prefix]
  
  cat(sprintf("\nProcessing Group %d (Prefix: %s)\n", i, current_prefix))
  cat(sprintf("  Samples in this group: %s\n", paste(samples_in_group, collapse=", ")))
  
  # このグループのBAMファイルを検索
  group_bams <- c()
  for (sample_id in samples_in_group) {
    # 完全一致に近いパターンマッチング - サンプルIDの前に境界または/があり、後ろに_R1などが続く
    pattern <- paste0("(^|/|[^A-Za-z0-9\\.])(", gsub("\\.", "\\\\.", sample_id), ")(_|$)")
    
    # デバッグ出力
    cat(sprintf("  Looking for pattern: %s\n", pattern))
    
    # 未使用のBAMファイルからマッチするものを検索
    available_bams <- setdiff(inputFiles, used_bams)
    matching_indices <- grep(pattern, available_bams)
    matching_bams <- available_bams[matching_indices]
    
    if (length(matching_bams) > 0) {
      cat(sprintf("  Found %d BAM file(s) for sample %s:\n", length(matching_bams), sample_id))
      for (bam in matching_bams) {
        cat(sprintf("    %s\n", basename(bam)))
      }
      group_bams <- c(group_bams, matching_bams)
      used_bams <- c(used_bams, matching_bams)  # 使用済みとしてマーク
    } else {
      cat(sprintf("  WARNING: No BAM files found for sample %s\n", sample_id))
    }
  }
  
  # 空でも全てのグループを保持
  group_files[[i]] <- group_bams
  
  # ファイルが見つからなかった場合は警告
  if (length(group_bams) == 0) {
    warning(sprintf("WARNING: No BAM files found for group %d (%s)", i, current_prefix))
  }
}

# 重複ファイルがないか確認
total_files <- sum(sapply(group_files, length))
if (total_files != length(unique(unlist(group_files)))) {
  warning("WARNING: Some BAM files appear to be assigned to multiple groups. Check the sample ID patterns in metadata.")
}

# group_sizesを計算
group_sizes <- sapply(group_files, length)

# 空のグループがあるか確認
empty_groups <- which(group_sizes == 0)
if (length(empty_groups) > 0) {
  warning(sprintf("WARNING: The following groups have no BAM files: %s", 
                 paste(paste0(empty_groups, " (", prefixes[empty_groups], ")"), collapse=", ")))
}

# 少なくとも1つのグループにBAMファイルがあるか確認
if (all(group_sizes == 0)) {
  stop("ERROR: No BAM files found for any group. Please check your metadata and BAM file names.")
}

# 空のグループがある場合、実行可能なグループだけを保持
if (length(empty_groups) > 0) {
  # ユーザーに空のグループを除外することを通知
  cat("\nRemoving empty groups from processing...\n")
  
  # 有効なグループのみを保持
  valid_groups <- which(group_sizes > 0)
  group_files <- group_files[valid_groups]
  prefixes <- prefixes[valid_groups]
  group_count <- length(valid_groups)
  group_sizes <- group_sizes[valid_groups]
  
  cat(sprintf("Continuing with %d valid groups: %s\n", 
             group_count, paste(prefixes, collapse=", ")))
}

# 全BAMファイルをフラット化
inputFiles <- unlist(group_files)

cat("\nSummary of groups based on metadata:\n")
cat(sprintf("Total groups: %d\n", group_count))
cat(sprintf("Group sizes: %s\n", paste(group_sizes, collapse=", ")))
cat(sprintf("Prefixes: %s\n", paste(prefixes, collapse=", ")))

### Generate sample labels
sampleLabels <- unlist(lapply(1:group_count, function(i) {
  paste0(prefixes[i], sprintf("%02d", 1:length(group_files[[i]])))
}))

### Generate mergeIndex
mergeIndex <- unlist(lapply(1:group_count, function(i) {
  rep(i, length(group_files[[i]]))
}))

# organism_name変数の設定
# パラメータから直接取得（優先）、もしくはシードファイルから抽出
if (!is.null(parsed_options$organismName)) {
  organism_name <- parsed_options$organismName
} else {
  # シードファイルから抽出した値を使用（バックアップ）
  organism_name <- sub("^BSgenome\\.([^\\.]+).*$", "\\1", package_name)
}

### TSSr object generation
myTSSr <- new("TSSr", genomeName = genome_name
                 ,inputFiles = inputFiles
                 ,inputFilesType = inputFilesType
                 ,sampleLabels = sampleLabels
                 ,sampleLabelsMerged = prefixes
                 ,mergeIndex = mergeIndex
                 ,refSource = refSource
                 ,organismName = organism_name)


# TSS calling
## options
# sequencingQualityThreshold: minimum sequencing quality score for a base to be considered (default:10, specified in integer)
# mappingQualityThreshold: minimum mapping quality score for a read to be considered (default:20, specified in integer)
# softClippingAllowed: (default: false), Note: If you ran STAR with the default alignment parameters (without using --alignEndsType Extend5pOfRead1), 
#                                              please ensure that you set softclippingAllowed = TRUE when running TSSr. 
##

getTSS(myTSSr, sequencingQualityThreshold = 20, mappingQualityThreshold = 20)


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


# Clustering TSSs to infer core promoters with "clusterTSS" function

# The "clusterTSS" function was designed to group neighboring TSSs into distinct TSS clusters (TCs), representing putative core promoters
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
           useMultiCore = TRUE, numCores = parsed_options$threads)


# Aggregating consensus TSS clusters with "consensusCluster" function.

# TSSr infers a set of consensus core promoters using the "consensusCluster" function to assign the same ID for TCs belong to the same core promoter
## options
# dis: (default: 50) Minimum distance between two peaks to be aggregated together into the same consensus cluster.
# useMultiCore: Logical indicating whether multiple cores are used or not. Default is FALSE.
# numCores: Number of cores. Default is NULL
##

consensusCluster(myTSSr, dis = 100, useMultiCore = TRUE, numCores = parsed_options$threads)


# export the consensus TCs to txt file
#
## options
# data: "tagClusters", "consensusClusters", "assigned", "unassigned". Default is "assigned"
##


# Assigning TCs to downstream genes

# Assigning TCs to downstream genes as their core promoters is required for annotation of the 5' boundaries of genomic features. 
# This process is also a prerequisite for further interrogations of regulated transcription initiation at the gene level. 
# TSSr offers the "annotateCluster" function to assign TCs to their downstream genes. The assignment of a TC to a gene is 
# based on the distance between the position of the dominant TSS of a TC and the annotated 5'ends of coding sequences (CDS) or transcripts.
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


# Export the annotated TCs in each group to txt files

# exported files prefix is each "sampleLabels"
# e.g. "sampleLabels".assignedClusters.txt

exportClustersTable(myTSSr, data = "assigned")
exportClustersTable(myTSSr, data = "consensusClusters")
exportClustersToBed(myTSSr, data = "consensusClusters")