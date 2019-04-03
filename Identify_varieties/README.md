# Variety confirmation 

This pipeline is used to identify the wheat cultivar in [field pathogenomics RNA-seq samples][https://academic.oup.com/gbe/article/9/12/3282/4644453]. It relies on the SNP makers developed with the [Breeders' 35K Axiom® array][ http://www.cerealsdb.uk.net/cerealgenomics/CerealsDB/axiom_download.php], which contains 35,143 SNPs selected to be informative across a wide range of hexaploid wheat accessions.

![FP-varietyconf][FP-varietyconf.png]

1. Copy the Identify_varieties directory

2. The reference contigs Reference_contigs/contigs_with_markers.fasta

3. Filter and trim the fastq reads if it's not already done

4. Align the reads against the reference_contigs using top_hat

   ``tophat -r 200 -o <sample>/top_hat_axiom contigs_with_marker $fastqfileR1 $fastqfileR2``

5. Sort and index the BAM files

   ``samtools sort <sample>/top_hat_axiom/accepted_hits.bam <sample>/BAM_files_axiom/accepted_hits_sorted``
   ``samtools index <sample>/BAM_files_axiom/accepted_hits_sorted.bam``

6. SNP calling

   ``samtools mpileup -f contigs_with_marker.fasta <sample>/BAM_files_axiom/accepted_hits_sorted.bam | gzip -9 -c > <sample>/accepted_hits_sorted.pileup.gz``

   ``gunzip -c <sample>/accepted_hits_sorted.pileup.gz | perl SCRIPTS/compsnps_pipe1_sampileup.py > <sample>/SNPs_axiom/<sample>\_SNP_ratios.txt``

   ``gunzip -c <sample>/accepted_hits_sorted.pileup.gz | java -jar VarScan.v2.3.9.jar mpileup2snp <sample>/accepted_hits_sorted.pileup --output-vcf 1 > <sample>/<sample>.vcf``

7. Obtain the reference tab file

   ``perl SCRIPTS/extract_2x_same_ref.pl <sample>/SNPs_axiom/<sample>\_SNP_ratios.txt > <sample>/<sample>\_Reference_greater_2x.tab``

8. To confirm the wheat variety, we need to run: 

   ``ruby Identify_varieties/SNP_markers.rb $sample $path``

``ruby Identify_varieties/decode_scores.rb $m ​$path``

where `$sample` is the name of the library and `$path` is the path where the library can be found

The file final_scores_`sample`_library.csv  will be generated.

Run the following R script to obtain the top score that corresponds to the identified variety: 

``LIB <- read.csv("final_scores_<sample>.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)
varieties <- colSums(LIB[,-1])
df <- data.frame(varieties)
write.csv(df, "<sample>.csv")
colMax(df)
varieties``

