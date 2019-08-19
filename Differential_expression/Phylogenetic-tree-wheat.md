# Field pathogenomics - Phylogenetic tree for wheat 

1. **Align using STAR to refseqv1.0 (whole genome)**

​    ``STAR --runThreadN 12 --runMode alignReads --genomeDir iwgsc_refseqv1.0_allchr.fa --readFilesIn $sample\_clean_R1.fastq $sample\_clean_R2.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $line/$line --sjdbGTFfile iwgsc_refseqv1.0_HighConf_2017Mar13.gff3 --sjdbGTFtagExonParentTranscript Parent``




wc -m < ${files[0]} >> $OUTPUT 

\# Write all the sequences to the final.phy file to generate a sequential PHYLIP file 

**for filename in $consensus/\*.filtered; do echo -n "$(basename $filename | cut -d '_' -f 1) "; cat $filename; echo ""; done >> $OUTPUT**



1. **Run tree**