
library(reshape2)
library(plyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggrepel)

d <- read.csv("~/Documents/PhD/Clusters-GO-annotation-FP-jan19/TCE-clust-early-upregulation/TCE-early-MF-goterms_by_cluster-level2-function.csv", header=TRUE, row.names=1)


####Heatmap using ggplot2#####


###### TCE selected genes clusters Heatmaps.  ######

d <- read.csv("~/Postdoc/clust/revigo-TCE-earlydownregulation-1314F22-OA-SA-SO-BP.csv", header=TRUE, row.names=1)

d_sorted<-d[order(-rowSums(d)),, drop = FALSE]
d.m <- melt(as.matrix(d_sorted))
head(d.m)

###rescale values
d.m <- ddply(d.m,.(Var2), transform,rescale = rescale(value, to=c(0, 100)))

#### convert 0s to NA 
d.m[d.m == 0] <- NA

pdf(paste("~/Postdoc/clust/revigo-TCE-earlyupregulation-1314F22-SA-OA-SO-BP-removeempty", ".pdf", sep="") ,onefile=TRUE, width = 9, height = 10)
p<-ggplot(d.m, aes(Var2, Var1)) + geom_tile(aes(fill = rescale), colour="grey") + scale_fill_distiller(palette = "Spectral", na.value="grey96") + theme(axis.text.x = element_text(angle = 90)) 
base_size <- 9
p + ggtitle("-log(pvalue) > 2)") + labs(x = "Cluster",y = "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0))
dev.off()

####  GO term enrichment barcharts  ####

d <- read.csv("~/Postdoc/1dpi-DEG-COM/gProfiler_taestivum_07-04-2020_09-11-37__intersections.csv", header=TRUE)
d<- na.omit(d)

d <- d[d$negative_log10_of_adjusted_p_value>5 | d$negative_log10_of_adjusted_p_value>5,]
d <- d[d$adjusted_p_value__C0<0.001 | d$adjusted_p_value__C0<0.001,]

d_bp<-d[d$source=="BP",]
d_mf<-d[d$source=="MF",]
d_cc<-d[d$source=="CC",]

barchart_bp <- ggplot(d_bp, aes(x = reorder(term, -log(adjusted_p_value__C0)), y =-log(adjusted_p_value__C0))) + geom_bar(stat = "identity",  width=.8, fill="red") +
  xlab("") + ylab("-log(adj p-value)") +theme_bw() + ggtitle("") + coord_flip()
 
barchart_cc <- ggplot(d_cc, aes(x = reorder(term, negative_log10_of_adjusted_p_value), y =negative_log10_of_adjusted_p_value)) + geom_bar(stat = "identity",  width=.8, fill="green") +
  xlab("") + ylab("-log(adj p-value)") +theme_bw() + ggtitle("")+ coord_flip()
 
barchart_mf <- ggplot(d_mf, aes(x = reorder(term, negative_log10_of_adjusted_p_value), y =negative_log10_of_adjusted_p_value)) + geom_bar(stat = "identity",  width=.8, fill="blue") +
  xlab("") + ylab("-log(adj p-value)") +theme_bw() + ggtitle("") + coord_flip()

filename="1dpi-1314-notF22-barchart2"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=20, height=8)
barchart_bp + barchart_cc + barchart_mf 
dev.off()


######  GO term enrichment Bubble plot  ####
d <- read.csv("~/Postdoc/1dpi-DEG-COM/1dpi-F22not1314-gProfiler_taestivum_forplot.csv", header=TRUE)

bubble <- ggplot(d, aes(x = term, y =negative_log10_of_adjusted_p_value, colour=source, size=intersection_size), label = d$term) + geom_point(aes(colour = factor(source), alpha=0.5)) +
  xlab("") + ylab("-log10 (Padj)") +theme_bw() + ggtitle("") + scale_size(range = c(1, 10)) + geom_text_repel(data=subset(d, negative_log10_of_adjusted_p_value > 4),aes(term, negative_log10_of_adjusted_p_value,label=stringr::str_wrap(term, 2)), check_overlap = TRUE, size=3)
bubble <- bubble + theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(0, 60))


filename="1dpi-1314-vs-F22-bubbles-2"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=12, height=12)
bubble
dev.off()

#Bubble plot function for enriched GO terms

d<-read.csv("Goterms-C0-1314OA-tab.txt")

bubble <- function(d, title){
  d <- d[d$adjusted_p_value<1 | d$adjusted_p_value<1,]
  bubble <- ggplot(d, aes(x = term, y =-log(adjusted_p_value), colour=source, size=intersection_size), label = d$term) + geom_point(aes(colour = factor(source), alpha=0.5)) +
    xlab("") + ylab("-log10 (Padj)") +theme_bw() + scale_size(range = c(1, 10)) + geom_text_repel(data=subset(d, -log(adjusted_p_value) > 20),aes(term, -log(adjusted_p_value),label=stringr::str_wrap(term, 2)), size=3)
  bubble1 <- bubble + ggtitle(title) + theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
}


C2 <- read.csv("Goterms-C2-1314OA-tab.txt", header=TRUE)
plotC2 <- bubble(C2, "C2")




#filename<-paste(gene, sep="")
filename <- "try"
pdf(paste(filename, ".pdf", sep="") ,onefile=TRUE, width=10, height=30)
plotC0 / plotC1 / plotC2 / plotC3 / plotC4 / plotC5 / plotC6 / plotC7 / plotC8 / plotC9 / plotC10 / plotC11 / plotC12 / plotC13
dev.off()
plotC1

library(ggforce)


install.packages("ggforce")         


library(ggplot2)
barchart <- function(d, tittle, value) {
  d <- d[d$adjusted_p_value<value | d$adjusted_p_value<value,]
barchart <- ggplot(d, aes(x = reorder(term, -log(adjusted_p_value)), y =-log(adjusted_p_value), fill=source)) + geom_bar(stat = "identity",  width=.8) +
  xlab("") + ylab("-log(adj p-value)") +theme_bw() + ggtitle(value) + coord_flip()
barchart
} 





