
library(ggplot2)
library(reshape2)
library(plyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(grid)

all_genes_df<-read.csv("~/Postdoc/TCE-kallisto-tpm-dataframe-updatedMay.csv", header=TRUE)

##calculate average for three replicates of each condition 
f1314.con.oak <- rowMeans(all_genes_df[,2:4])
f22.con.oak <- rowMeans(all_genes_df[,2:4])
f22.con.sa <- rowMeans(all_genes_df[,5:7])
f1314.con.sa <- rowMeans(all_genes_df[,5:7])
f22.con.so <- rowMeans(all_genes_df[,8:10])
f1314.con.so <- rowMeans(all_genes_df[,8:10])

f22.oak.t1 <-rowMeans(all_genes_df[,11:13])
f22.sa.t1  <-rowMeans(all_genes_df[,14:16])
f22.so.t1 <-rowMeans(all_genes_df[,17:19])

f22.oak.t3 <-rowMeans(all_genes_df[,20:22])
f22.sa.t3 <-rowMeans(all_genes_df[,23:25])
f22.so.t3  <-rowMeans(all_genes_df[,26:28])

f22.oak.t7 <-rowMeans(all_genes_df[,29:31])
f22.sa.t7 <-rowMeans(all_genes_df[,32:34])
f22.so.t7 <-rowMeans(all_genes_df[,35:37])

f22.oak.t11  <-rowMeans(all_genes_df[,38:40])
f22.sa.t11 <-rowMeans(all_genes_df[,41:43])
f22.so.t11 <-rowMeans(all_genes_df[,44:46])

f1314.oak.t1 <-rowMeans(all_genes_df[,47:49])
f1314.sa.t1  <-rowMeans(all_genes_df[,50:53])
f1314.so.t1 <-rowMeans(all_genes_df[,54:56])

f1314.oak.t3 <-rowMeans(all_genes_df[,57:59])
f1314.sa.t3 <-rowMeans(all_genes_df[,59:61])
f1314.so.t3  <-rowMeans(all_genes_df[,61:63])

f1314.oak.t7 <-rowMeans(all_genes_df[,64:66])
f1314.sa.t7 <-rowMeans(all_genes_df[,67:69])
f1314.so.t7 <-rowMeans(all_genes_df[,70:72])

f1314.oak.t11  <-rowMeans(all_genes_df[,73:75])
f1314.sa.t11 <-rowMeans(all_genes_df[,76:78])
f1314.so.t11 <-rowMeans(all_genes_df[,79:81])

f22.sa <- data.frame(all_genes_df$target_id, "control"=f22.con.sa,  "t1"=f22.sa.t1, "t3"=f22.sa.t3, "t7"=f22.sa.t7,"t11"=f22.sa.t11)
f22.sol <- data.frame(all_genes_df$target_id, "control"=f22.con.so, "t1"=f22.so.t1, "t3"=f22.so.t3, "t7"=f22.so.t7,"t11"=f22.so.t11)
f22.oak <- data.frame(all_genes_df$target_id, "control"=f22.con.oak, "t1"=f22.oak.t1, "t3"=f22.oak.t3, "t7"=f22.oak.t7,"t11"=f22.oak.t11)

f1314sa <- data.frame(all_genes_df$target_id, "control"=f1314.con.sa , "t1"=f1314.sa.t1, "t3"=f1314.sa.t3, "t7"=f1314.sa.t7,"t11"=f1314.sa.t11)
f1314sol <- data.frame(all_genes_df$target_id, "control"=f1314.con.so, "t1"=f1314.so.t1, "t3"=f1314.so.t3, "t7"=f1314.so.t7,"t11"=f1314.so.t11)
f1314oak <- data.frame(all_genes_df$target_id, "control"=f1314.con.oak, "t1"=f1314.oak.t1, "t3"=f1314.oak.t3, "t7"=f1314.oak.t7,"t11"=f1314.oak.t11)

###PLOT FOR GENE LIST####

#####Compare F22-Santiago and 13/14-Santiago in same plot 

genes <- read.table("~/Postdoc/TCE-paper/Results-DEG-vsCON/NLR-specific-toSAF22.txt")

for (gene in genes$V1){
  dfsant<-f22.sa[all_genes_df$target_id == gene,]
  dfsant[["isolate"]] <- c("F22")
  dfsant1314<-f1314sa[all_genes_df$target_id == gene,]
  dfsant1314[["isolate"]] <- c("1314")
  newdff <- rbind(dfsant, dfsant1314)
  newdff <- melt(newdff, id.vars=c("all_genes_df.target_id", "isolate"))
  g <- ggplot(newdff, aes(y=value,x=as.numeric(variable), colour=isolate)) +geom_line(size=2) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle("")+ xlim("Control","1","3","7","11") 
  filename<-paste(gene, sep="")
  filename
  pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE)
  grid.arrange(g, top=textGrob(gene, gp=gpar(fontsize=15,font=8, face="bold")))
  dev.off()
}


#####For all cultivar-isolate pairs, F22 separated from 13/14

genes <- read.table("~/Postdoc/TCE-paper/selected-genes-chloroplast.txt")

for (gene in genes$V1){
  dfoak1314<-f1314oak[all_genes_df$target_id == gene,]
  dfoak1314[["variety"]] <- c("Oakley")
  dfsols1314<-f1314sol[all_genes_df$target_id == gene,]
  dfsols1314[["variety"]] <- c("Solstice")
  dfsant1314<-f1314sa[all_genes_df$target_id == gene,]
  dfsant1314[["variety"]] <- c("Santiago")
  df1314 <- rbind(dfoak1314, dfsols1314, dfsant1314)
  dfm1314melted <- melt(df1314, id.vars=c("all_genes_df.target_id", "variety"))
  dfoak<-f22.oak[all_genes_df$target_id == gene,]
  dfoak[["variety"]] <- c("Oakley")
  dfsols<-f22.sol[all_genes_df$target_id == gene,]
  dfsols[["variety"]] <- c("Solstice")
  dfsant<-f22.sa[all_genes_df$target_id == gene,]
  dfsant[["variety"]] <- c("Santiago")
  dff22 <- rbind(dfoak, dfsols, dfsant)
  dfmf22melted <- melt(dff22, id.vars=c("all_genes_df.target_id", "variety"))
  p <- ggplot(dfm1314melted, aes(y=value,x=as.numeric(variable),colour=variety)) +geom_line(size=2) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + scale_color_manual(values=c('orangered2','lightskyblue1', 'green')) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text=element_blank()) + ggtitle("") + xlim("Control","1","3","7","11")
  g <- ggplot(dfmf22melted, aes(y=value,x=as.numeric(variable), colour=variety)) +geom_line(size=2) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text=element_blank()) + ggtitle("")+ xlim("Control","1","3","7","11") 
  filename<-paste(gene, sep="")
  filename
  pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE)
  grid.arrange(g,p, top=textGrob(gene, gp=gpar(fontsize=15,font=8, face="bold")))
  dev.off()
}


###FOR UNIQUE GENE ID####

gene="TraesCS2A02G247300.2"

##build plot for F22 
dfoak<-f22.oak[all_genes_df$target_id == gene,]
dfoak[["variety"]] <- c("Oakley")
dfsols<-f22.sol[all_genes_df$target_id == gene,]
dfsols[["variety"]] <- c("Solstice")
dfsant<-f22.sa[all_genes_df$target_id == gene,]
dfsant[["variety"]] <- c("Santiago")
dff22 <- rbind(dfoak, dfsols, dfsant)

dfmf22melted <- melt(dff22, id.vars=c("all_genes_df.target_id", "variety"))

g <- ggplot(dfmf22melted, aes(y=value,x=as.numeric(variable), colour=variety)) +geom_line(size=2) + theme_bw() +
  ylab("tpm") + xlab("Days post-inoculation") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle("")+ xlim("Control","1","3","7","11") 
g
##build plot for 13/14
dfoak1314<-f1314oak[all_genes_df$target_id == gene,]
dfoak1314[["variety"]] <- c("Oakley")
dfsols1314<-f1314sol[all_genes_df$target_id == gene,]
dfsols1314[["variety"]] <- c("Solstice")
dfsant1314<-f1314sa[all_genes_df$target_id == gene,]
dfsant1314[["variety"]] <- c("Santiago")
df1314 <- rbind(dfoak1314, dfsols1314, dfsant1314)
dfm1314melted <- melt(df1314, id.vars=c("all_genes_df.target_id", "variety"))

p <- ggplot(dfm1314melted, aes(y=value,x=as.numeric(variable),colour=variety)) +geom_line(size=2) + theme_bw() +
  ylab("tpm") + xlab("Days post-inoculation") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle("") + xlim("Control","1","3","7","11")
p

### make output file 
filename<-paste(gene, sep="")
filename
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE)
grid.arrange(g,p, top=textGrob(gene, gp=gpar(fontsize=15,font=8, face="bold")))
dev.off()


### Box plot 1 dpi comparison F22 vs 13/14
data<-read.csv("sig-2-F22specific-GO0072593ROS-foldchange.csv", header=TRUE)
mtdf<-melt(data)
write.csv(mtdf, "sig-2-F22specific-GO0072593ROS-foldchange.csv")
datamelted<-read.csv("sig-2-F22specific-GO0072593ROS-foldchange.csv", header=TRUE)

p <- ggplot(data = datamelted, aes(x=variable, y=value)) + ylim(0,6) +ylab("Fold-change(tpm +1)") + geom_boxplot(aes(fill=variable))+ theme_bw() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p + facet_wrap( ~ gene, scales="free")
