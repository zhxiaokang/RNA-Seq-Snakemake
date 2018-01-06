library(openxlsx)

dataBaP.low <- read.xlsx('table.xlsx', sheet = 1, colNames = T, rowNames = T)
dataBaP.high <- read.xlsx('table.xlsx', sheet = 2, colNames = T, rowNames = T)
dataEE2.low <- read.xlsx('table.xlsx', sheet = 3, colNames = T, rowNames = T)
dataEE2.high <- read.xlsx('table.xlsx', sheet = 4, colNames = T, rowNames = T)
dataMix.low <- read.xlsx('table.xlsx', sheet = 5, colNames = T, rowNames = T)
dataMix.high <- read.xlsx('table.xlsx', sheet = 6, colNames = T, rowNames = T)

barplot(t(as.matrix(dataBaP.low)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))

barplot(t(as.matrix(dataBaP.high)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))

barplot(t(as.matrix(dataEE2.low)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))

barplot(t(as.matrix(dataEE2.high)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))

barplot(t(as.matrix(dataMix.low)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))

barplot(t(as.matrix(dataMix.high)),ylim=c(0,1),xlab='P-Value',ylab="Proportion",main="Proportions",col=rainbow(3),font=2,cex.lab=1.4)
legend("topleft",legend=c("DESeq2_uniq","Overlap","edgeR_uniq"),pch=15,col=rainbow(3))
