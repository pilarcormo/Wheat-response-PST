# Field pathogenomics - Phylogenetic tree for wheat 

All scripts and full pipelines to make a phylogenetic tree using field RNA-seq samples can be be found in [https://github.com/vbuens/Field_Pathogenomics/blob/master/Tree_pipeline.sh](https://github.com/vbuens/Field_Pathogenomics/blob/master/Tree_pipeline.sh). These are the adaptations made for wheat: 

### 1. Filter using bbduk

	fastqfileR1=$line/$line\_R1.fastq;
	fastqfileR2=$line/$line\_R2.fastq;
    gunzip -c $line/*$line*\_R1.fastq.gz > $line/$fastqfileR1;
    gunzip -c $line/*$line*\_R2.fastq.gz > $line/$fastqfileR2;
    /nbi/software/testing/bbtools/37.68/x86_64/bbduk.sh -Xmx1g t=12 in1=$line/$fastqfileR1 in2=$line/$fastqfileR2 out1=$line/$line\_clean_R1.fastq out2=$line/$line\_clean_R2.fastq minlen=36 qtrim=rl trimq=10 ktrim=r k=25 mink=11 ref=$path/adapters/TruSeq3-PE-2.fa hdist=1

### 2. Align using STAR to refseqv1.0 (whole genome)

	STAR --runThreadN 12 --runMode alignReads --genomeDir $genome --readFilesIn $line/$line\_clean_R1.fastq $line/$line\_clean_R2.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $line/$line --sjdbGTFfile $genome/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3 --sjdbGTFtagExonParentTranscript Parent

### 3. Consensus files

	samtools mpileup -uf $genome/iwgsc_refseqv1-dna.fa $home/$line/*.out.bam | bcftools call -c | $scripts/vcfutils.pl vcf2fq > $line/$line\_cns.fastq
	seqtk seq -a  $line/$line\_cns.fastq > $line/$line\_consensus.fasta
	rm -rf $line/$line\_cns.fastq
	tr '[:lower:]' '[:upper:]' < $line/$line\_consensus.fasta > $line/$line\_consensus_uc.fasta

### 4. Concatenation and filtering 

	FILES="*consensus_uc.fasta"
	shopt -s nullglob
	for f in $FILES;
	do
  		python $home/SCRIPTS/phylogenetic_trees/sort_fasta.py -f $f > $consensus/$(basename $f).sorted
	done 
	# With all the files ordered by gene and placed in the sorted folder, filter out which genes are relevant (have enough information)
	##### codon_from_fasta has multiple filtering parameters, the most relevant are -l (minimum percentage of known bases in a sequence for acceptance) -s (minimum number of accepted samples percentage)
	python $home/SCRIPTS/phylogenetic_trees/codon_from_fasta.py -d $consensus -l 20 -s 20 -c codon123

	# Select all the gene filtered sequences
	pattern="$consensus/*.filtered"
	files=( $pattern )

	# Write the PHYLIP header to the final.phy file
	echo -n ${#files[@]} > $OUTPUT
	echo -n " " >> $OUTPUT

	wc -m < ${files[0]} >> $OUTPUT
	# Write all the sequences to the final.phy file to generate a sequential PHYLIP file which can be used to run the tre. 
	for filename in $consensus/*.filtered; do echo -n "$(basename $filename | cut -d '_' -f 1) "; cat $filename; echo ""; done >> $OUTPUT
