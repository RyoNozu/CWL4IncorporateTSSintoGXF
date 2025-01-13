# Import required files & TSS calling
inputFiles <- list.files(path = "./", pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
inputFilesType <- "bamPairedEnd" # set “inputFilesType” as “bamPairedEnd” for paired-end BAM files, and as "TSStable" if the input file is a TSS table 

args <- commandArgs(trailingOnly=TRUE)
refSource <- args[1] # コマンドライン引数からrefSourceを取得
prefix1 <- args[2] # コマンドライン引数から1つ目のprefixを取得
prefix2 <- args[3] # コマンドライン引数から2つ目のprefixを取得

# グループごとにファイルを分ける
group1_files <- inputFiles[1:(length(inputFiles)/2)]
group2_files <- inputFiles[(length(inputFiles)/2 + 1):length(inputFiles)]

# サンプルラベルを生成
sampleLabels_group1 <- paste0(prefix1, sprintf("%02d", 1:length(group1_files)))
sampleLabels_group2 <- paste0(prefix2, sprintf("%02d", 1:length(group2_files)))
sampleLabels <- c(sampleLabels_group1, sampleLabels_group2)

myTSSr <- new("TSSr", genomeName = "BSgenome.Htrimaculatus.inhouse.Htriv1"
                 ,inputFiles = inputFiles
                 ,inputFilesType = inputFilesType
                 ,sampleLabels = sampleLabels
                 ,sampleLabelsMerged = c("group1","group2")
                 ,mergeIndex = c(rep(1, length(group1_files)), rep(2, length(group2_files)))
                 ,refSource = refSource
                 ,organismName = "Halichoeres trimaculatus")

myTSSr

getTSS(myTSSr, sequencingQualityThreshold=20, mappingQualityThreshold=20)
