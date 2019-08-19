library(DESeq2)
library(ggplot2)
library("RColorBrewer")
library("genefilter")
library("gplots")

directory <- "~/Differential_expression/"

#directory <- "~/Differential_expression/Htseq_files/"

#cond<-c("Oakley","Cordiale","Oakley","Cordiale","Santiago","Santiago","Oakley","Santiago","Santiago","Oakley","Oakley","Oakley","Santiago","Oakley","Oakley","None","Oakley","Cordiale","Cordiale","Santiago")
#type<-c("Oakley","Cordiale","Oakley","Cordiale","Santiago","Santiago","Oakley","Santiago","Santiago","Oakley","Oakley","Oakley","Santiago","Oakley","Oakley","None","Oakley","Cordiale","Cordiale","Santiago")

cond<-c("Oakley","Cordiale","Oakley","Cordiale","Santiago","Santiago","Oakley")
type<-c("Oakley","Cordiale","Oakley","Cordiale","Santiago","Santiago","Oakley")


file<-paste(directory, "htseq_tab.txt", sep="") #####Only 2015 samples 


#Prepare the table
countsTable<-read.table(file,header=T)
rownames(countsTable)<-countsTable$gene
countsTable<-countsTable[,-1] 
rlog = rlogTransformation


#############
#COMPARE ALL#
#############
localTable<-countsTable
localCond<-cond

colData<-data.frame(condition=factor(localCond), type=factor(type))
localCond
head(colData)
#The levels need to be set
dds<-DESeqDataSetFromMatrix(countData=localTable,colData,formula(~condition))

ddsMF <- DESeq(dds)

res <- results(ddsMF)

res_filtered <- subset(resOrdered , res$pvalue < 0.05)
write.csv(res_filtered, "deseq_output.csv")
resOrdered <- res[order(res$padj),]
rld <- rlog(dds)
rld
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
#pca = prcomp(t(assay(x)))

plotMA(dds,ylim=c(-5,5),main=name)

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) + theme_bw() +
  geom_point(size=7) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



##########HEATMAP COMPARING ALL
sampleDists <- dist(t(assay(rld)))  
sampleDistMatrix <- as.matrix(sampleDists)

library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)


library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ), 50 )


all_genes <- assay(rld)[topVarGenes,]
all_genes_df <- data.frame(all_genes)

assay(rld)[topVarGenes,]

means_13 <- rowMeans(all_genes_df[,3:7])
means_24 <- rowMeans(all_genes_df[,2:4])
means_56 <- rowMeans(all_genes_df[,5:6]) 


mean_dataframe <- data.frame("Oakley"=means_13, "Cordiale"=means_24, "Santiago"=means_56)


library("gplots")

heatmap.2(as.matrix(mean_dataframe), scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins=c(10,25),lhei = c(0.15, 0.95))


###############HEATMAP ROWSUMS

#check distribution of row sums
quantile(rowSums(countsTable))
#create a workable set
data_subset <- countsTable[rowSums(countsTable)>50000,]
nrow(data_subset)
data_matrix <- data.matrix(data_subset)
heatmap(data_matrix)

heatmap.2(data_matrix,col=brewer.pal(11,"RdBu"),scale="row", trace="none", margins=c(10,25))



########################
##PAIRWISE COMPARISON##
#######################

Cordiale<-c("LIB18683", "LIB19042")
Oakley<-c("LIB18260", "LIB19021", "LIB19530")
Santiago<-c("LIB19340", "LIB19332")

#####Oakley vs Santiago
name<-paste("Oakley", "vs", "Santiago", sep=" ")
filename<-paste("Oakley", "vs", "Santiago",sep="")

localTable<-countsTable[,c(Oakley, Santiago)]
localCond<-c(rep("Oakley", 3), rep("Santiago", 2))
colData<-data.frame(condition=factor(localCond))
localCond
head(colData)
#The levels need to be set
dds<-DESeqDataSetFromMatrix(countData=localTable,colData,formula(~condition))
dds<-DESeq(dds)
res <- results(dds)
res_filtered <- subset(res, res$padj < 0.05)
resOrdered <- res_filtered[order(res_filtered$padj),]
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
pdf(paste("~/Differential_expression","/", filename, "MA.pdf", sep="") ,onefile=TRUE)
plotMA(dds,ylim=c(-5,5),main=name)

###write output csv
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file=paste("~/Differential_expression","/" ,filename, ".csv", sep=""))
#####

topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ), 50)

heatmap.2(assay(rld)[topVarGenes,], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), margins=c(10,25))


#####


####Santiago vs Cordiale
name<-paste("Santiago", "vs", "Cordiale", sep=" ")
filename<-paste("Santiago", "vs", "Cordiale",sep="")

localTable<-countsTable[,c(Santiago, Cordiale)]
localCond<-c(rep("Santiago", 2), rep("Cordiale", 2))
colData<-data.frame(condition=factor(localCond))
localCond
head(colData)
#The levels need to be set
dds<-DESeqDataSetFromMatrix(countData=localTable,colData,formula(~condition))
dds<-DESeq(dds)

res <- results(dds)
res_filtered <- subset(res, res$padj < 0.05)
resOrdered <- res_filtered[order(res_filtered$padj),]

rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
pdf(paste("~/Differential_expression","/", filename, "MA.pdf", sep="") ,onefile=TRUE)
plotMA(dds,ylim=c(-5,5),main=name)

###write output csv
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file=paste("~/Differential_expression","/" ,filename, ".csv", sep=""))
#####

topVarGenes2 <- head(order(rowVars(assay(rld)), decreasing=TRUE ), 100)

heatmap.2(assay(rld)[topVarGenes2,], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins=c(10,25))



#####

#####Oakley vs Cordiale
name<-paste("Oakley", "vs", "Cordiale", sep=" ")
filename<-paste("Oakley", "vs", "Cordiale",sep="")

localTable<-countsTable[,c(Oakley, Cordiale)]
localCond<-c(rep("Oakley", 3), rep("Cordiale", 2))
colData<-data.frame(condition=factor(localCond))
localCond
head(colData)
#The levels need to be set
dds<-DESeqDataSetFromMatrix(countData=localTable,colData,formula(~condition))
dds<-DESeq(dds)

res <- results(dds)
res_filtered <- subset(res , res$padj < 0.05)
resOrdered <- res_filtered[order(res_filtered$padj),]
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ), 100)
library("gplots")

heatmap.2(assay(rld)[topVarGenes,], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins=c(10,25))


pdf(paste("~/Differential_expression","/", filename, "MA.pdf", sep="") ,onefile=TRUE)
plotMA(dds,ylim=c(-5,5),main=name)

resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file=paste("~/Differential_expression","/" ,filename, ".csv", sep=""))

#####

