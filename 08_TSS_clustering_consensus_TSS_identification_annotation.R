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

clusterTSS(myTSSr, method = "peakclu", peakDistance=100, extensionDistance=20, 
           localThreshold = 0.02, clusterThreshold = 1, 
           useMultiCore=TRUE, numCores=10)


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
                upstream=5000, upstreamOverlap = 500, downstream = 0)

### Export the annotated TCs in each group to txt files
# exported files prefix is each "sampleLabels" (refer to 07_ R script)
# e.g. "sampleLabels".assignedClusters.txt

exportClustersTable(myTSSr, data = "assigned")