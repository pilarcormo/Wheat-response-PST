
library(tximportData)
library(tximport)
library(DESeq2)
library("RColorBrewer")
library(ggplot2)

####Set working directory where files are
setwd("/Users/morenop/Postdoc/Kallisto-results")
directory <- "/Users/morenop/Postdoc/Kallisto-results"

libs1dpi <- grep("kallisto-", list.files(directory),value=T)
kal_dirs <- sapply(libs1dpi, function(id) file.path(directory, id, "abundance.h5"))
all(file.exists(kal_dirs))

file<-paste("htseq-TCE.txt", sep="")
cond<-read.table(file, header=TRUE, row.names=1)

####tximport kallisto files into counts 
txi.kallisto <- tximport(kal_dirs, type = "kallisto", txOut = TRUE)

####tximport counts to DESeq2 
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto.1dpi,cond,design = ~isolate+cultivar)

dds1dpi <- DESeq(ddsTxi)
rld <- vst(dds1dpi)

##PCA
pcaData <- plotPCA(rld, intgroup=c("isolate", "cultivar"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste("PCA", ".pdf", sep="") ,onefile=TRUE, width=10, height=10)
ggplot(pcaData, aes(PC1, PC2, color=timepoint, shape=cultivar)) + theme_bw() +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


