
library(ggplot2)
library(reshape2)
library(plyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyverse)

all_genes_df<-read.csv("TCE-kallisto-tpm-dataframe-updatedMay.csv", header=TRUE)

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
f1314sol <- data.frame(all_genes_df$target_id,"control"=f1314.con.so, "t1"=f1314.so.t1, "t3"=f1314.so.t3, "t7"=f1314.so.t7,"t11"=f1314.so.t11)
f1314oak <- data.frame(all_genes_df$target_id, "control"=f1314.con.oak, "t1"=f1314.oak.t1, "t3"=f1314.oak.t3, "t7"=f1314.oak.t7,"t11"=f1314.oak.t11)


###dataframe for subset already
f22.sa <- f22.sa %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
f22.sol <- f22.sol %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
f22.oak <- f22.oak %>% filter_at(vars(starts_with("control")), any_vars(. > 1))

f1314sa <- f1314sa %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
f1314sol <- f1314sol %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
f1314oak <- f1314oak %>% filter_at(vars(starts_with("control")), any_vars(. > 1))

dfmf22oak <- melt(f22.oak)
dfmf22oak <- ddply(dfmf22oak,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sol <- melt(f22.sol)
dfmf22sol <- ddply(dfmf22sol,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sant <- melt(f22.sa)
dfmf22sant <- ddply(dfmf22sant,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314oak <- melt(f1314oak)
dfm1314oak <- ddply(dfm1314oak,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf1314sol <- melt(f1314sol)
dfmf1314sol <- ddply(dfmf1314sol,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314sant <- melt(f1314sa)
dfm1314sant <- ddply(dfm1314sant,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))


####input gene ID list #####
genes <- read.table("~/Postdoc/TCE-paper/Immunity/GO0000302mart_export.txt")

newdfoak <- df[FALSE,]
newdfsol <- df[FALSE,]
newsant <- df[FALSE,]
newdfoak1314 <- df[FALSE,]
newdfsol1314 <- df[FALSE,]
newsant1314 <- df[FALSE,]

#Generate expression dataset for gene list 
for (gene in genes$V1){
  dfoak<-f22.oak[f22.oak$all_genes_df.target_id == gene,]
  newdfoak <- rbind(newdfoak, dfoak)
  dfsols<-f22.sol[f22.sol$all_genes_df.target_id == gene,]
  newdfsol <- rbind(newdfsol, dfsols)
  dfsant<-f22.sa[f22.sa$all_genes_df.target_id == gene,]
  newsant <- rbind(newsant, dfsant)
  dfoak1314<-f1314oak[f1314oak$all_genes_df.target_id == gene,]
  newdfoak1314 <- rbind(newdfoak1314, dfoak1314)
  dfsols1314<-f1314sol[f1314sol$all_genes_df.target_id == gene,]
  newdfsol1314 <- rbind(newdfsol1314, dfsols1314)
  dfsant1314<-f1314sa[f1314sa$all_genes_df.target_id == gene,]
  newsant1314 <- rbind(newsant1314, dfsant1314)
}

newsant <- newsant %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
newdfsol <- newdfsol %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
newdfoak <- newdfoak %>% filter_at(vars(starts_with("control")), any_vars(. > 1))

newsant1314 <- newsant1314 %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
newdfsol1314 <- newdfsol1314 %>% filter_at(vars(starts_with("control")), any_vars(. > 1))
newdfoak1314 <- newdfoak1314 %>% filter_at(vars(starts_with("control")), any_vars(. > 1))

dfmf22oak <- melt(newdfoak)
dfmf22oak <- ddply(dfmf22oak,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sol <- melt(newdfsol)
dfmf22sol <- ddply(dfmf22sol,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf22sant <- melt(newsant)
dfmf22sant <- ddply(dfmf22sant,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314oak <- melt(newdfoak1314)
dfm1314oak <- ddply(dfm1314oak,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfmf1314sol <- melt(newdfsol1314)
dfmf1314sol <- ddply(dfmf1314sol,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
dfm1314sant <- melt(newsant1314)
dfm1314sant <- ddply(dfm1314sant,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))
################



gene_plot<-function(df, titlename, colour){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=all_genes_df.target_id)) +geom_line(size=0,alpha=0) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun = "median", geom = "line", color = colour, size = 1.2)
  g <- g + xlim("Control","1","3","7","11")
  g
}
gene_plot1314<-function(df, titlename, colour){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=all_genes_df.target_id)) +geom_line(size=0, alpha=0) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun = "median", geom = "line", color = colour, size = 1.2)
  g <- g + xlim("Control","1","3","7","11")
  g
}

oak<-gene_plot(unique(dfmf22oak), "Oakley F22", "indianred3")
sol<-gene_plot(unique(dfmf22sol), "Solstice F22", "deepskyblue4")
san<-gene_plot(unique(dfmf22sant), "Santiago F22", "darkolivegreen4")
oak1314<-gene_plot1314(unique(dfm1314oak), "Oakley 13/14","indianred1")
sol1314<-gene_plot1314(unique(dfmf1314sol), "Solstice 13/14", "deepskyblue1")
san1314<-gene_plot1314(unique(dfm1314sant), "Santiago 13/14", "darkolivegreen3")

setwd("~/Postdoc/Specific-go-terms-plots/")
filename="GO0050789-allgenes-allcultivars-filter1tpm"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=10, height=10)
oak + sol + san + oak1314 + sol1314 + san1314
dev.off()

san1314GO0015979<-gene_plot1314(unique(dfm1314sant), "Susceptible", "indianred3")
sanGO0015979<-gene_plot(unique(dfmf22sant), "Resistant", "darkolivegreen4")
sanGO0015979 + san1314GO0015979

filename="GO0050789-allgenes-Santiago-filter1tpm"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=10, height=5)
sanGO0015979 + san1314GO0015979
dev.off()
