library(ggplot2)
library(reshape2)
library(plyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyverse)

setwd("~/Postdoc/TCE-paper/")
all_genes_df<-read.csv("~/Postdoc/TCE-paper/TCE-kallisto-tpm-dataframe-updatedMay.csv", header=TRUE)

##control
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

f22.sa <- data.frame("f22.sacontrol"=f22.con.sa,  "f22.sat1"=f22.sa.t1, "f22.sat3"=f22.sa.t3, "f22.sat7"=f22.sa.t7,"f22.sat11"=f22.sa.t11)
f22.sol <- data.frame("f22.solcontrol"=f22.con.so, "f22.solt1"=f22.so.t1, "f22.solt3"=f22.so.t3, "f22.solt7"=f22.so.t7,"f22.solt11"=f22.so.t11)
f22.oak <- data.frame("f22.oakcontrol"=f22.con.oak, "f22.oakt1"=f22.oak.t1, "f22.oakt3"=f22.oak.t3, "f22.oakt7"=f22.oak.t7,"f22.oakt11"=f22.oak.t11)

f1314sa <- data.frame("f1314sacontrol"=f1314.con.sa , "f1314sat1"=f1314.sa.t1, "f1314sat3"=f1314.sa.t3, "f1314sat7"=f1314.sa.t7,"f1314sat11"=f1314.sa.t11)
f1314sol <- data.frame("f1314solcontrol"=f1314.con.so, "f1314solt1"=f1314.so.t1, "f1314solt3"=f1314.so.t3, "f1314solt7"=f1314.so.t7,"f1314solt11"=f1314.so.t11)
f1314oak <- data.frame("f1314oakcontrol"=f1314.con.oak, "f1314oakt1"=f1314.oak.t1, "f1314oakt3"=f1314.oak.t3, "f1314oakt7"=f1314.oak.t7,"f1314oakt11"=f1314.oak.t11)

merged <- data.frame(all_genes_df$target_id,f22.sa, f22.sol, f22.oak, f1314sa, f1314sol, f1314oak) 

####input gene ID list and output filenames#####
genes <- read.table("~/Postdoc/TCE-paper/Immunity/GO0045087mart_export.txt")
filename1="scaled-GO0045087-ROS-allgenes-allcultivars-filter1tpm"
filename2="scaled-GO0045087-ROS-1314-specificgenesonly25-allgenes-Santiago-filter1tpm"

newmergedsubset <- df[FALSE,]

#Generate expression dataset for gene list 

for (gene in genes$V1){
  mergedsubset <- merged[merged$all_genes_df.target_id == gene,]
  newmergedsubset <- rbind(newmergedsubset, mergedsubset)
}

newmergedsubset <- newmergedsubset %>% filter_at(vars(starts_with("f22.sacontrol")), any_vars(. > 1))
meltednewmergedsubset <- melt(newmergedsubset)
meltednewmergedsubset <- ddply(meltednewmergedsubset,.(all_genes_df.target_id), transform, rescale = rescale(value, to=c(0, 1)))


dfmf22oak <- dplyr::filter(meltednewmergedsubset, !grepl("sol|sa|f1314",variable))
dfmf22sol <- dplyr::filter(meltednewmergedsubset, !grepl("oak|sa|f1314",variable))
dfmf22sant <- dplyr::filter(meltednewmergedsubset, !grepl("sol|oak|f1314",variable))
dfm1314oak <- dplyr::filter(meltednewmergedsubset, !grepl("sol|sa|f22",variable))
dfmf1314sol <- dplyr::filter(meltednewmergedsubset, !grepl("oak|sa|f22",variable))
dfm1314sant <- dplyr::filter(meltednewmergedsubset, !grepl("sol|oak|f22",variable))


gene_plot<-function(df, titlename, colour){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=all_genes_df.target_id)) +geom_line(size=0,alpha=0) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun = "median", geom = "line", color = colour, size = 1.2)
  g <- g + ylim(0,1)
  g
}

gene_plot1314<-function(df, titlename, colour){
  g <- ggplot(df, aes(y=rescale,x=as.numeric(variable), colour=all_genes_df.target_id.x)) +geom_line(size=0,alpha=0) + theme_bw() +
    ylab("tpm") + xlab("Days post-inoculation") + theme(legend.position="none") + scale_color_grey() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + ggtitle(titlename) + stat_summary(fun = "median", geom = "line", color = colour, size = 1.2)
  g <- g + ylim(0,1)+ xlim(6,10)
}

oak<-gene_plot(unique(dfmf22oak), "Oakley F22", "indianred3")
sol<-gene_plot(unique(dfmf22sol), "Solstice F22", "deepskyblue4")
san<-gene_plot(unique(dfmf22sant), "Santiago F22", "darkolivegreen4")
oak1314<-gene_plot(unique(dfm1314oak), "Oakley 13/14","indianred1")
sol1314<-gene_plot(unique(dfmf1314sol), "Solstice 13/14", "deepskyblue1")
san1314<-gene_plot(unique(dfm1314sant), "Santiago 13/14", "darkolivegreen3")

pdf(paste(filename1, ".pdf", sep="") ,onefile=TRUE, width=10, height=10)
oak + sol + san + oak1314 + sol1314 + san1314
dev.off()

f22<-gene_plot(unique(dfmf22sant), "Santiago F22-Resistant", "darkolivegreen4")
f1314<-gene_plot(unique(dfm1314sant), "Santiago 1314-Susceptible", "indianred3")
f22 + f1314

pdf(paste(filename2, ".pdf", sep="") ,onefile=TRUE, width=10, height=5)
f22 + f1314
dev.off()
