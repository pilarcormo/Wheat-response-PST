
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(patchwork)

###Function#########################
cluster_plot<-function(df, titlename, line){
  g <- ggplot(df, aes(y=value,x=as.numeric(variable), colour=Genes)) +geom_line(size=0.1) + theme_bw() +
    ylab("tpm") + xlab("dpi") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10))  + stat_summary(fun.y = "median", geom = "line", color = line, size = 1.2)
  g <- g + xlim("Control","1","3","7","11") + ggtitle(titlename)
  g
} 
##################################

####Example code 
##input files are the output from running Clust and can be retrieved for each specific cluster by running: 
## awk -F '\t' '{print $1}' Clusters_Objects.tsv > C0.txt 
## grep -f C0.txt Input_files_and_params/Data/Data_processed.tsv >> C0-datapoints.txt 


cluster <- read.table("C0-datapoints.txt", header=TRUE)
C0 <- melt(cluster)

cluster <- read.table("C1-datapoints.txt", header=TRUE)
C1 <- melt(cluster)

cluster <- read.table("C2-datapoints.txt", header=TRUE)
C2 <- melt(cluster)

cluster <- read.table("C3-datapoints.txt", header=TRUE)
C3 <- melt(cluster)

cluster <- read.table("C4-datapoints.txt", header=TRUE)
C4 <- melt(cluster)

cluster <- read.table("C5-datapoints.txt", header=TRUE)
C5 <- melt(cluster)

cluster <- read.table("C6-datapoints.txt", header=TRUE)
C6 <- melt(cluster)

cluster <- read.table("C7-datapoints.txt", header=TRUE)
C7 <- melt(cluster)

cluster <- read.table("C8-datapoints.txt", header=TRUE)
C8 <- melt(cluster)

cluster <- read.table("C9-datapoints.txt", header=TRUE)
C9 <- melt(cluster)

cluster <- read.table("C10-datapoints.txt", header=TRUE)
C10 <- melt(cluster)

cluster <- read.table("C11-datapoints.txt", header=TRUE)
C11 <- melt(cluster)

cluster <- read.table("C12-datapoints.txt", header=TRUE)
C12 <- melt(cluster)

cluster <- read.table("C13-datapoints.txt", header=TRUE)
C13 <- melt(cluster)


plotC0 <- cluster_plot(C0, "Cluster I (2069 genes)", "blue")
plotC1 <- cluster_plot(C1, "Cluster IX (127 genes)", "black")
plotC2 <- cluster_plot(C2, "Cluster V (82 genes)", "darkgreen")
plotC3 <- cluster_plot(C3, "Cluster IV (213 genes)", "red")
plotC4 <- cluster_plot(C4, "Cluster II (1419 genes)", "red")
plotC5 <- cluster_plot(C5, "Cluster III (194 genes)", "red")
plotC6 <- cluster_plot(C6, "Cluster X (216 genes)", "black")
plotC7 <- cluster_plot(C7, "Cluster VI (529 genes)", "darkgreen")
plotC8 <- cluster_plot(C8, "Cluster VII (183 genes)", "darkgreen")
plotC9 <- cluster_plot(C9, "Cluster VIII (139 genes)", "darkgreen")
plotC10 <- cluster_plot(C10, "Cluster XI (63 genes)", "black")


filename<-paste(gene, sep="")
filename <- "cluster-colours"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=10, height=10)
plotC0 + plotC4 + plotC5 + plotC3 + plotC2 + plotC7 + plotC8 + plotC9  + plotC1  + plotC6  + plotC10
dev.off()

##