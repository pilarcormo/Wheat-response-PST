# Field pathogenomics - host differential gene expression

This pipeline is used to find differentially expressed genes in rust-infected wheat using [field pathogenomics RNA-seq samples](https://academic.oup.com/gbe/article/9/12/3282/4644453). It's adapted for the [recently released](http://science.sciencemag.org/content/361/6403/eaar7191/tab-figures-data) wheat annotation [IWGSC Refseq v1.1](https://wheat-urgi.versailles.inra.fr/Seq-Repository/Annotations). 

![FP-DEG](/Users/morenop/Documents/Wheat-response-PST/Differential_expression/FP-kallisto.png)





1. Read quality and adapter trimming 

   ```
   bbduk.sh -Xmx1g t=12 in1=$fastqfileR1 in2=$fastqfileR2 out1=$sample\_clean_R1.fastq out2=$sample\_clean_R2.fastq minlen=36 qtrim=rl trimq=10 ktrim=r k=25 mink=11 ref=adapters/TruSeq3-PE-2.fa hdist=1
   ```

2. Calculate transcript abundaces using [Kallisto](https://pachterlab.github.io/kallisto/)

   ```bash
   kallisto index -i IWGSC_v1.1_HLC.idx ​IWGSC_v1.1_HLC.fasta``
   
   kallisto quant -i WGSC_v1.1_HLC.idx -o kallisto_iwgsc_​$sample  -b 100 $  $sample\_clean_R1.fastq $sample\_clean_R2.fastq
   ```

   

3. Run [tximport](https://github.com/mikelove/tximport) in R to get count data from Kallisto output files 

   ```R
   library(tximportData)
   library(tximport)
   
   ### Find Kallisto output files ###
   
   libs <- grep("kallisto_iwgsc", list.files(directory),value=T)
   kal_dirs <- sapply(libs, function(id) file.path(directory, id, "abundance.h5"))
   all(file.exists(kal_dirs))
   
   #### tximport kallisto files into counts ###
   
   txi.kallisto <- tximport(kal_dirs, type = "kallisto", txOut = TRUE)
   counts<-txi.kallisto$counts
   ncol(counts)
   
   ## Make sure all count values are integers and not numeric format ###
   sapply(counts, class)
   
   ### Write csv output table with counts ##
   write.csv(counts, "table-counts-fp18.csv")
   ```

   4. Import metadata file with replicates information

      ```R
      file<-paste(directory, "../replicates.txt", sep="")
      cond<-read.table(file, header=TRUE, row.names=1)
      
      ### Check that metadata file contains the libraries in the same order as counts table 
      
      all(rownames(cond) %in% colnames(counts))
      all(rownames(cond) == colnames(counts))
      counts <- counts[, rownames(cond)]
      all(rownames(cond) == colnames(counts))
      rownames(cond)
      ```

      5. Run RUVseq

   ``` R
   #### re-order data and make sure the values are integers ###
   
   filter <- apply(countsTableOrder, 1, function(x) length(x[x>5])>=2)
   filtered <- countsTableOrder[filter,]
   filtered[,-1] <- sapply(filtered[,-1], as.integer)
   
   ### Combine count data with metadata to generate set ###
   set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(cond$Variety, row.names=colnames(filtered)))
   
   ### Make RLE and PCA plots before normalisation ###        
   condition<-cond$Variety
   colors <- brewer.pal(5, "Set2")
   plotRLE(set, outline=FALSE, ylim=c(-3, 3), col=colors[condition],las=2, legend = TRUE)
   plotPCA(set, col=colors[cond$Variety], cex=1,label=FALSE) 
   cond$Variety
   
   ### Set 2 - RUVg ###               
   design <-model.matrix(~cond$Variety, data=pData(set))
   y <- DGEList(counts=counts(set), group=cond$Variety)
   y <- calcNormFactors(y, method="upperquartile")
   y <- estimateGLMCommonDisp(y, design)
   y <- estimateGLMTagwiseDisp(y, design)
   fit <- glmFit(y, design)
   lrt <- glmLRT(fit, coef=2)
   top <- topTags(lrt, n=nrow(set))$table
   empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
   set2 <- RUVg(set, empirical, k=1)
   pData(set2)
   plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[cond],las=2)
   plotPCA(set2, col=colors[cond$Variety], cex=1)
                   
   ### Set 3 - RUVs method to estimate the factors of unwanted variation using replicate/negative control samples for which the covariates of interest are constant.###  
                   
   #First, we need to construct a matrix specifying the replicates.
   differences <- makeGroups(cond$Variety)
   differences
   genes <- rownames(filtered)
   set3 <- RUVs(set, genes, k=3, differences)
   data<- pData(set3)
   plotRLE(set3, outline=FALSE, ylim=c(-3, 3), col=colors[condition],las=2, legend = TRUE)
   colors <- brewer.pal(3, "Set2")
   plotPCA(set3, col=colors[cond$Variety], label=TRUE, cex=1)
                   
   ### Set 4 - RUVr ###
   design <- model.matrix(~condition, data=pData(set))
   y <- DGEList(counts=counts(set), group=condition)
   y <- calcNormFactors(y, method="upperquartile")
   y <- estimateGLMCommonDisp(y, design)
   y <- estimateGLMTagwiseDisp(y, design)
   fit <- glmFit(y, design)
   res <- residuals(fit, type="deviance")
   set4 <- RUVr(set, genes, k=1, res)
   plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[condition],las=2)
   plotPCA(set4, col=colors[condition], cex=1, las=2)              
                          
   ```

   

6. Import normalised counts from RUVseq to DESeq2 



Clust](https://github.com/BaselAbujamous/clust)





https://github.com/tanghaibao/goatools





TILLING lines with high impact mutations on gene of interest

[impact_mut.rb](https://github.com/pilarcormo/Wheat-response-PST/blob/master/TILLING/impact_mut.rb)