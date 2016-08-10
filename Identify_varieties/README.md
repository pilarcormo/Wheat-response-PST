###Variety confirmation

- Copy the ```Identify_varieties``` directory

- The reference contigs are ```Reference_contigs/contigs_with_markers.fasta```

- Filter and trim the fastq reads if it's not already done. 
- Align the reads to the **reference_contigs** using top_hat
- Sort and index the BAM files
- mpileup (bam to pileup)

```
m=$1  ##name of the library
samtools mpileup -f $work/$Genome\.fasta $m/accepted_hits_sorted.bam | gzip -9 -c > $m/accepted_hits_sorted.pileup.gz
```
- SNP calling using VarScan (Yoy need to have VarScan.v2.3.9.jar)

```
source jre-6.45
line=$1 ##name of the library
	cd $line
	gunzip -c accepted_hits_sorted.pileup.gz > accepted_hits_sorted.pileup
	java -jar VarScan.v2.3.9.jar mpileup2snp accepted_hits_sorted.pileup --output-vcf 1 > $line.vcf 
	cd ..
```

- Obtain the reference tab file

```
m=$1  ##name of the library
mkdir $m/SNPs
gunzip -c $m/accepted_hits_sorted.pileup.gz | perl $work/SCRIPTS/compsnps_pipe1_sampileup.py > $m/SNPs/$m\_SNP_ratios.txt
perl $work/SCRIPTS/extract_2x_same_ref.pl $m/SNPs/$m\_SNP_ratios.txt > $m/$m\_Reference_greater_2x.tab
```

- Confirm de variety:

```
ruby Identify_varieties/SNP_markers_with_vcf.rb $m 
```
From this script you will obtain the file ``snp_markers_<library>.csv`` that contains 4 fields (name of the markers, contig, position, nt in the sample, score: 2 for homozygous, 1 for heterozygous and 0 for reference). 

```
ruby Identify_varieties/decode_scores.rb $m 
```
This script decodes the previous score asigned to each position for each marker. From this script you will obtain the file ```final_scores_<library>.csv``` that contains all the varieties and all the markers and a score of 0, 0.5 or 1 for each marker and variety. 

Need to specify as arguments: 

- ```$m``` Name of the library 



###R
1- With the output cvs file ```final_scores_<library>.csv``` we need to add together all the values in the columns to calculate the score for each variety:

```
LIB <- read.csv("final_scores_<library>.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)
varieties <- colSums(LIB[,-1])
df <- data.frame(varieties)
colMax(df)
varieties
```
We will see the maximum score and all the varieties and their score assigned. Look for the variety that matches the highest number and that will be the predicted variety for the sample. 

2- To make the plot, run the following after step 1: 

```
csv <- write.csv(df, "try.csv", quote=FALSE)
for_plot <- read.csv("~/try.csv")
df2 <- data.frame(for_plot$X, for_plot$varieties)
n<-dim(df)[1]
df2<-df2[1:(n-1),]
caracol <- ggplot(df2, aes(x = reorder(for_plot.X, for_plot.varieties), y = for_plot.varieties)) + geom_bar(width = 1, stat = "identity") +coord_polar() +theme_bw() + ylab(" ") + xlab(" ") + theme(plot.title = element_text(size = 30), axis.text=element_text(size=30), axis.title=element_text(size=50,face="bold")) + theme(axis.text.x = element_text(angle = 360/(2*pi)*rev( pi/2 + seq( pi/21, 2*pi-pi/21, len=21))+ 360/(2*pi)*c( rep(0, 10),rep(pi,10), rep(0,10))))
#------------------------
##optional to save
#------------------------
ggsave(name, device = "pdf")
dev.off()
```