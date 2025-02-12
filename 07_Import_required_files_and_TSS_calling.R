# Import required files & TSS calling
inputFiles <- list.files(path = "./", pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
inputFilesType <- "bamPairedEnd"  # set “inputFilesType” as “bamPairedEnd” for paired-end BAM files, and as "TSStable" if the input file is a TSS table 

args <- commandArgs(trailingOnly = TRUE)
refSource <- args[1] # コマンドライン引数からrefSourceを取得
group_count <- as.integer(args[2]) # コマンドライン引数からグループ数を取得
group_sizes <- as.integer(args[3:(2 + group_count)]) # コマンドライン引数から各グループのファイル数を取得
prefixes <- args[(3 + group_count):(2 + 2 * group_count)] # コマンドライン引数から各グループのprefixを取得

# グループごとにファイルを分ける
start_index <- 1
group_files <- list()
for (i in 1:group_count) {
  end_index <- start_index + group_sizes[i] - 1
  group_files[[i]] <- inputFiles[start_index:end_index]
  start_index <- end_index + 1
}

# サンプルラベルを生成
sampleLabels <- unlist(lapply(1:group_count, function(i) {
  paste0(prefixes[i], sprintf("%02d", 1:length(group_files[[i]])))
}))

# mergeIndexを生成
mergeIndex <- unlist(lapply(1:group_count, function(i) {
  rep(i, length(group_files[[i]]))
}))

myTSSr <- new("TSSr", genomeName = "BSgenome.Htrimaculatus.inhouse.Htriv1"
                 ,inputFiles = inputFiles
                 ,inputFilesType = inputFilesType
                 ,sampleLabels = sampleLabels
                 ,sampleLabelsMerged = paste0("group", 1:group_count)
                 ,mergeIndex = mergeIndex
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
