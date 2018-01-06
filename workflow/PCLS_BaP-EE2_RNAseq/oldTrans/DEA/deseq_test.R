library(DESeq2)

# load the count file
# args <- commandArgs(TRUE)
# countTable <- read.table(args[1], header = T)


countTable <- read.table('./Mix_high_VS_control', header = T)
rownames(countTable) <- countTable[,1]
countData <- countTable[,-1] # remove Gene_id column to get pure count table

# save gene list in geneList for extracting gene names later
geneList <- rownames(countData)

# get the sample id
samples <- colnames(countData)

# get the fish group
numSample <- length(samples)
groupTreat <- seq(1, numSample/2)  # index different fish group
groupControl <- seq(1, numSample/2)  # make sure that treat and control group are in pairs
fishGroup <- as.character(c(groupTreat, groupControl))

# generate the experiment group
expGroup <- c(rep('treat', numSample/2), rep('control', numSample/2))

# get the colData mixing the above 3 data
colData <- data.frame(samples, fishGroup, expGroup)

# creat the DESeqDataSet
dds_paird <- DESeqDataSetFromMatrix(countData, colData, design = ~ fishGroup + expGroup)

# perform DEA
dea_paired <- DESeq(dds_paird, fitType = 'parametric', test = 'LRT', reduced = 'T')

# result analysis
res_paired <- results(dea_paired)

# p-adjusted value
padj_paired <- res_paired$padj

# significantly expressed genes
sig_paired_local <- geneList[which(padj_paired < 0.05)]

# write.csv(res_paired, args[2])