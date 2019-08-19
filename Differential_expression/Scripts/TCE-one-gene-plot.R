
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

##Example gene 
gene="TraesCS4A01G346900.1"


colors <- brewer.pal(3, "Set2")

###Find gene of interest in dataset for F22-infected samples
dfoak<-f22.oak[f22.oak$countsTable.gene == gene,]
dfoak[["variety"]] <- c("Oakley")
dfsols<-f22.sol[f22.sol$countsTable.gene == gene,]
dfsols[["variety"]] <- c("Solstice")
dfsant<-f22.sa[f22.sa$countsTable.gene == gene,]
dfsant[["variety"]] <- c("Santiago")
bindf22 <- rbind(dfoak, dfsols, dfsant)
dfmf22 <- melt(bindf22, id.vars=c("countsTable.gene", "variety"))

###Find gene of interest in dataset for 13/14-infected samples
dfoak1314<-f1314oak[f1314oak$countsTable1314.gene == gene,]
dfoak1314[["variety"]] <- c("Oakley")
dfsols1314<-f1314sol[f1314sol$countsTable1314.gene == gene,]
dfsols1314[["variety"]] <- c("Solstice")
dfsant1314<-f1314sa[f1314sa$countsTable1314.gene == gene,]
dfsant1314[["variety"]] <- c("Santiago")
bindf1314 <- rbind(dfoak1314, dfsols1314, dfsant1314)
dfm1314 <- melt(bindf1314, id.vars=c("countsTable1314.gene", "variety"))

#make expression plots 
plotf22<- ggplot(dfmf22, aes(y=value,x=as.numeric(variable), colour=variety)) +geom_line(size=2) + theme_bw() +
  ylab("tpm") + xlab("Days post-inoculation") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle("")+ xlim("Control","1","3","7","11") 
plot1314 <- ggplot(dfm1314, aes(y=value,x=as.numeric(variable),colour=variety)) +geom_line(size=2) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle("") + xlim("Control","1","3","7","11")

#save expression plots 
filename<-paste(gene, sep="")
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE)
grid.arrange(plotf22,plot1314, top=textGrob(gene, gp=gpar(fontsize=15,font=8, face="bold")))
dev.off()

