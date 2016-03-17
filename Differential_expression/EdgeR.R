
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)

x <- read.delim("~/Differential_expression/htseq_tab.txt",row.names="gene")
group <- factor(c(1, 2, 1, 2, 3, 3, 1)) #### 1 = oakley, 2 = cordiale, 3 = santiago


y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)

design <- model.matrix(~group)

y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)


qlf_oak_cord <- glmQLFTest(fit,coef=2) ###COMPARING OAKLEY VS CORDIALE
qlf_cord_san <- glmQLFTest(fit,contrast=c(0,-1,1)) ####COMPARING CORDIALE VS SANTIAGO
qlf_oak_san <- glmQLFTest(fit,coef=3) ###COMPARING OAKLEY VS SANTIAGO


plotMDS(y)


top_oak_san<-topTags(qlf_oak_san, n=23537, sort.by = "PValue") ### n= number of genes 
top_cord_san<-topTags(qlf_cord_san, n=23537, sort.by = "PValue") ### n= number of genes 
top_oak_cord<-topTags(qlf_oak_cord, n=23537, sort.by = "PValue") ### n= number of genes 

df_oak_san<-data.frame(top_oak_san)
df_cord_san <- data.frame(top_cord_san)
df_oak_cord <-data.frame(top_oak_cord)

###filter the data to have only significant genes 

de_genes_oak_cord <- subset(df_oak_cord, df_oak_cord$FDR < 0.05)
de_genes_cord_san <- subset(df_cord_san, df_cord_san$FDR < 0.05)
de_genes_oak_san <- subset(df_oak_san, df_oak_san$FDR < 0.05)

###write output
write.csv(de_genes_oak_san, file ="output_edgeR_oakvssan.csv")
write.csv(de_genes_cord_san, file ="output_edgeR_cordvssan.csv")
write.csv(de_genes_oak_cord, file ="output_edgeR_oakvscord.csv")




