library(edgeR)
library(statmod)

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

# Put the data into a DGEList object
y <- DGEList(counts = countTable[,2:(length(samples)+1)], genes = countTable[,1])

# Normalization

# compute the library sizes
y$samples$lib.size<-colSums(y$counts)

# define the group
#y$samples$group <- factor(factor(c(rep('treat', numSample/2), rep('control', numSample/2))))

# use Gene_id as row names
rownames(y$counts) <- rownames(y$genes) <- y$genes$genes
y$genes$genes <- NULL

# TMM normalization is applied to this dataset to account for compositional difference between the libraries
y <- calcNormFactors(y)

# The design matrix

# get the fish group
numSample <- length(samples)
groupTreat <- seq(1, numSample/2)  # index different fish group
groupControl <- seq(1, numSample/2)  # make sure that treat and control group are in pairs
fishGroup <- factor(c(groupTreat, groupControl))

# generate the experiment group
expGroup <- factor(c(rep('treat', numSample/2), rep('control', numSample/2)))

design <- model.matrix(~fishGroup+expGroup)
rownames(design) <- colnames(y)

# Estimating the dispersion

# estimate the NB dispersion for the dataset
y <- estimateDisp(y, design, robust = TRUE)

# Differential expression

# determine differentially expressed genes
# fit genewise glms
fit <- glmFit(y, design)

# conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes
lrt <- glmLRT(fit)

# differentially expressed genes
de <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
toptag <- topTags(lrt, n = length(geneList), p.value = 0.05)
deg <- rownames(toptag$table)

