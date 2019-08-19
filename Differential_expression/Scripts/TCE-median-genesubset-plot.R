setwd("~/Documents/PhD/KALLISTO/")

library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)

###################################Create dataset for F22-infected samples ###################################
countsTable<-read.csv("F22.csv", header=TRUE)
all_genes_df <- data.frame(countsTable)

f22.con.sa <- rowMeans(all_genes_df[,5:7])
f22.con.so <- rowMeans(all_genes_df[,8:9])
f22.con.oak <- rowMeans(all_genes_df[,2:4])
f22.sa.t1 <- rowMeans(all_genes_df[,13:15])
f22.so.t1 <- rowMeans(all_genes_df[,16:18])
f22.oak.t1 <-rowMeans(all_genes_df[,10:12])
f22.sa.t11 <-rowMeans(all_genes_df[,22:24])
f22.so.t11 <-rowMeans(all_genes_df[,25:27])
f22.oak.t11  <-rowMeans(all_genes_df[,19:21])
f22.sa.t3 <-rowMeans(all_genes_df[,31:33])
f22.so.t3 <-rowMeans(all_genes_df[,34:36])
f22.oak.t3 <-rowMeans(all_genes_df[,28:30])
f22.sa.t7 <-rowMeans(all_genes_df[,40:42])
f22.so.t7 <-rowMeans(all_genes_df[,43:44])
f22.oak.t7 <-rowMeans(all_genes_df[,37:39])


f22.sa <- data.frame(countsTable$gene, "control"=f22.con.sa,  "t1"=f22.sa.t1, "t3"=f22.sa.t3, "t7"=f22.sa.t7,"t11"=f22.sa.t11)
f22.sol <- data.frame(countsTable$gene, "control"=f22.con.so, "t1"=f22.so.t1, "t3"=f22.so.t3, "t7"=f22.so.t7,"t11"=f22.so.t11)
f22.oak <- data.frame(countsTable$gene, "control"=f22.con.oak, "t1"=f22.oak.t1, "t3"=f22.oak.t3, "t7"=f22.oak.t7,"t11"=f22.oak.t11)


###################################Create dataset for 13/14-infected samples ###################################

countsTable1314<-read.csv("1314_2.csv", header=TRUE)

all_genes_df <- data.frame(countsTable1314)

f1314.con.sa <- rowMeans(all_genes_df[,5:7])
f1314.con.so <- rowMeans(all_genes_df[,8:9])
f1314.con.oak <- rowMeans(all_genes_df[,2:4])
f1314.sa.t1 <- rowMeans(all_genes_df[,13:15])
f1314.so.t1 <- rowMeans(all_genes_df[,16:18])
f1314.oak.t1 <-rowMeans(all_genes_df[,10:12])
f1314.sa.t11 <-rowMeans(all_genes_df[,22:24])
f1314.so.t11 <-rowMeans(all_genes_df[,25:27])
f1314.oak.t11  <-rowMeans(all_genes_df[,19:21])
f1314.sa.t3 <-rowMeans(all_genes_df[,31:33])
f1314.so.t3 <-rowMeans(all_genes_df[,34:36])
f1314.oak.t3 <-rowMeans(all_genes_df[,28:30])
f1314.sa.t7 <-rowMeans(all_genes_df[,40:42])
f1314.so.t7 <-rowMeans(all_genes_df[,43:44])
f1314.oak.t7 <-rowMeans(all_genes_df[,37:39])

f1314sa <- data.frame(countsTable1314$gene, "control"=f1314.con.sa , "t1"=f1314.sa.t1, "t3"=f1314.sa.t3, "t7"=f1314.sa.t7,"t11"=f1314.sa.t11)
f1314sol <- data.frame(countsTable1314$gene,"control"=f1314.con.so, "t1"=f1314.so.t1, "t3"=f1314.so.t3, "t7"=f1314.so.t7,"t11"=f1314.so.t11)
f1314oak <- data.frame(countsTable1314$gene, "control"=f1314.con.oak, "t1"=f1314.oak.t1, "t3"=f1314.oak.t3, "t7"=f1314.oak.t7,"t11"=f1314.oak.t11)


###Example gene list 
genes <- read.table("GO0023014-proteinphosphorylation.txt")

newdfoak <- df[FALSE,]
newdfsol <- df[FALSE,]
newsant <- df[FALSE,]
newdfoak1314 <- df[FALSE,]
newdfsol1314 <- df[FALSE,]
newsant1314 <- df[FALSE,]

#Generate expression dataset for gene list 
for (gene in genes$V1){
  dfoak<-f22.oak[f22.oak$countsTable.gene == gene,]
  newdfoak <- rbind(newdfoak, dfoak)
  dfsols<-f22.sol[f22.sol$countsTable.gene == gene,]
  newdfsol <- rbind(newdfsol, dfsols)
  dfsant<-f22.sa[f22.sa$countsTable.gene == gene,]
  newsant <- rbind(newsant, dfsant)
  dfoak1314<-f1314oak[f1314oak$countsTable1314.gene == gene,]
  newdfoak1314 <- rbind(newdfoak1314, dfoak1314)
  dfsols1314<-f1314sol[f1314sol$countsTable1314.gene == gene,]
  newdfsol1314 <- rbind(newdfsol1314, dfsols1314)
  dfsant1314<-f1314sa[f1314sa$countsTable1314.gene == gene,]
  newsant1314 <- rbind(newsant1314, dfsant1314)
}

###Normalise expression (0-1)
dfmf22oak <- melt(newdfoak)
dfmf22oak <- ddply(dfmf22oak,.(countsTable.gene), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sol <- melt(newdfsol)
dfmf22sol <- ddply(dfmf22sol,.(countsTable.gene), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sant <- melt(newsant)
dfmf22sant <- ddply(dfmf22sant,.(countsTable.gene), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314oak <- melt(newdfoak1314)
dfm1314oak <- ddply(dfm1314oak,.(countsTable1314.gene), transform, rescale = rescale(value, to=c(0, 1)))
dfmf1314sol <- melt(newdfsol1314)
dfmf1314sol <- ddply(dfmf1314sol,.(countsTable1314.gene), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314sant <- melt(newsant1314)
dfm1314sant <- ddply(dfm1314sant,.(countsTable1314.gene), transform, rescale = rescale(value, to=c(0, 1)))


#######functions to plot gene subset 
gene_plot<-function(df, titlename){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=countsTable.gene)) +geom_line(size=0.1) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun.y = "median", geom = "line", color = "black", size = 1.2)
  g <- g + xlim("Control","1","3","7","11")
  g
}
gene_plot1314<-function(df, titlename){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=countsTable1314.gene)) +geom_line(size=0.1) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun.y = "median", geom = "line", color = "black", size = 1.2)
  g <- g + xlim("Control","1","3","7","11")
  g
}

#######Plot gene subsets expression
oak<-gene_plot(unique(dfmf22oak), "Oakley")
sol<-gene_plot(unique(dfmf22sol), "Solstice")
san<-gene_plot(unique(dfmf22sant), "Santiago")
oak1314<-gene_plot1314(unique(dfm1314oak), "Oakley")
sol1314<-gene_plot1314(unique(dfmf1314sol), "Solstice")
san1314<-gene_plot1314(unique(dfm1314sant), "Santiago")

##combine them in a single plot 
one<-grid.arrange(oak,sol,san,top=textGrob("F22", gp=gpar(fontsize=15,font=8, face="bold")))
two<-grid.arrange(oak1314,sol1314,san1314,top=textGrob("13/14", gp=gpar(fontsize=15,font=8, face="bold")))

##Save plot image to pdf 
filename<-paste("uniq_photosynthesis-blast-refseq", sep="")
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE)
grid.arrange(one, two, ncol=2)
dev.off()
