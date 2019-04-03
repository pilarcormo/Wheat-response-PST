######CUFFDIFF
source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
install.packages("sqldf")
########

dbFile<-paste('~/Differential_expression/Oak-San-Cord_v1/cufflinks_2015_2',"/cuffData.db", sep="")


library(cummeRbund)
library(sqldf)

cuff<-readCufflinks('~/Differential_expression/Oak-San-Cord_v1/cufflinks_2015_2') ###create a database out of cuffdiff output files. 


gene_diff_data <- diffData(genes(cuff))


sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
diffGenes<-getGenes(cuff, sig_gene_data$gene_id)

sig_gene_data
diffGenes

names<-featureNames(diffGenes)


#myGeneIds<- c("Oakley","Santiago", "Cordiale")
#myGenes <- getGenes(cuff_data, myGeneIds)
#myGeneIds

####A volcano plot is a scatter plot that also identifies differentially expressed genes (by color) between samples

disp<-dispersionPlot(genes(cuff))
disp

genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
isoforms.scv

###distribution ofRNA-seq read counts (fpkm)
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep
brep<-csBoxplot(genes(cuff),replicates=T)
brep

##scatter okit 
scatter<-csScatterMatrix(genes(cuff))
scatter
scatter<-csScatterMatrix(genes(cuff),replicates=T)
scatter

###volcano plot
v<-csVolcanoMatrix(genes(cuff))
v
v<-csVolcanoMatrix(genes(cuff), replicates=T)
tab<-getSigTable(cuff)

##overview of significant features 
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05) 

###Distance matrix
myDistHeat<-csDistHeat(genes(cuff))
myDistHeat

#######Dimensionality reduction
###PCA
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2", replicates=T)
genes.PCA 

genes.MDS.rep<-MDSplot(genes(cuff),replicates=T) 
genes.MDS.rep

####Gene clusters 
ic <- csCluster(diffGenes,  k=10) 
head(ic)

csClusterPlot(clustering, pseudocount=1.0,logMode=TRUE,drawSummary=TRUE,sumFun=mean_cl_boot)

csHeatmap(diffGenes, cluste='both')
ic
featureNames(ic)

head(ic$cluster)  
icp<-csClusterPlot(ic) 
icp
write.table(ic$clustering, "~/Differential_expression/clusters_april.txt")

b<-expressionBarplot(myGenes) 
b 
myGenes.spec<-csSpecificity(diffGenes)  

##Heatmap
h<-csHeatmap(myGenes,cluster='both') 
h