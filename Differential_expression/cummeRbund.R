######CUFFDIFF

source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")


dbFile<-paste('~/Differential_expression/cuffdiff',"/cuffData.db", sep="")


library(cummeRbund)
library(sqldf)

cuff<-readCufflinks('~/Differential_expression/cuffdiff') ###create a database out of cuffdiff output files. 


gene_diff_data <- diffData(genes(cuff))


sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
diffGenes<-getGenes(cuff, sig_gene_data$gene_id)
names<-featureNames(diffGenes)

g <- genes(cuff)

csDensity(genes(cuff)) ###distribution ofRNA-seq read counts (fpkm)
csVolcanoMatrix(g) 
csScatter(g)


myGeneIds<- c("Oakley","Santiago", "Cordiale")
myGenes <- getGenes(cuff_data, myGeneIds)
myGeneIds

####A volcano plot is a scatter plot that also identifies differentially expressed genes (by color) between samples

disp<-dispersionPlot(genes(cuff))
disp

genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
isoforms.scv
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
brep<-csBoxplot(genes(cuff),replicates=T)
brep
scatter<-csScatterMatrix(genes(cuff))
scatter
 
scatter<-csScatterMatrix(genes(cuff),replicates=T)
scatter
v<-csVolcanoMatrix(genes(cuff))
v
v<-csVolcanoMatrix(genes(cuff), replicates=T)
v
tab<-getSigTable(cuff)

