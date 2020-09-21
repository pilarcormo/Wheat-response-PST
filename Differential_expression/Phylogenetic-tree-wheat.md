# Field pathogenomics - Phylogenetic tree for wheat 

1. Filter using bbduk

   fastqfileR1=$line/$line\_R1.fastq;
     fastqfileR2=$line/$line\_R2.fastq;
    gunzip -c $line/*$line*\_R1.fastq.gz > $line/$fastqfileR1;
    gunzip -c $line/*$line*\_R2.fastq.gz > $line/$fastqfileR2;
    /nbi/software/testing/bbtools/37.68/x86_64/bbduk.sh -Xmx1g t=12 in1=$line/$fastqfileR1 in2=$line/$fastqfileR2 out1=$line/$line\_clean_R1.fastq out2=$line/$line\_clean_R2.fastq minlen=36 qtrim=rl trimq=10 ktrim=r k=25 mink=11 ref=$path/adapters/TruSeq3-PE-2.fa hdist=1

2. Align using STAR to refseqv1.0 (whole genome)

    STAR --runThreadN 12 --runMode alignReads --genomeDir $genome --readFilesIn $line/$line\_clean_R1.fastq $line/$line\_clean_R2.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix $line/$line --sjdbGTFfile $genome/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3 --sjdbGTFtagExonParentTranscript Parent

3. Consensus files

samtools mpileup -uf $genome/iwgsc_refseqv1-dna.fa $home/$line/*.out.bam | bcftools call -c | $scripts/vcfutils.pl vcf2fq > $line/$line\_cns.fastq
seqtk seq -a  $line/$line\_cns.fastq > $line/$line\_consensus.fasta
rm -rf $line/$line\_cns.fastq
tr '[:lower:]' '[:upper:]' < $line/$line\_consensus.fasta > $line/$line\_consensus_uc.fasta

4. Concatenation and filtering 

FILES="*consensus_uc.fasta"
shopt -s nullglob
for f in $FILES;
do
  python $home/SCRIPTS/phylogenetic_trees/sort_fasta.py -f $f > $consensus/$(basename $f).sorted
done 
##### With all the files ordered by gene and placed in the sorted folder, filter out which genes are relevant (have enough information)
##### codon_from_fasta has multiple filtering parameters, the most relevant are -l (minimum percentage of known bases in a sequence for acceptance) -s (minimum number of accepted samples percentage)
python $home/SCRIPTS/phylogenetic_trees/codon_from_fasta.py -d $consensus -l 20 -s 20 -c codon123

##### Select all the gene filtered sequences
pattern="$consensus/*.filtered"
files=( $pattern )

##### Write the PHYLIP header to the final.phy file
echo -n ${#files[@]} > $OUTPUT
echo -n " " >> $OUTPUT

wc -m < ${files[0]} >> $OUTPUT
##### Write all the sequences to the final.phy file to generate a sequential PHYLIP file
for filename in $consensus/*.filtered; do echo -n "$(basename $filename | cut -d '_' -f 1) "; cat $filename; echo ""; done >> $OUTPUT

5. Run tree

/tgac/software/testing/raxml/8.0.20/x86_64/raxmlHPC-PTHREADS-SSE3 -T 10 -s $OUTPUT -m GTRGAMMA -n tree.newick -p 100

6. Run bootstraps for tree

name=('tree')
consensus=('boots')
mkdir $consensus
cd $consensus;
mkdir seeds
for step in $(seq $SLURM_ARRAY_TASK_ID) ; do
	seed=$RANDOM;
done

echo $seed > seeds/$SLURM_ARRAY_TASK_ID.seed.txt
/tgac/software/testing/raxml/8.0.20/x86_64/raxmlHPC-PTHREADS-SSE3 -T 10 -s ../alignment.fasta -m GTRGAMMA -n $name\_3rd_codon_tree_boots_$SLURM_ARRAY_TASK_ID.phy -p 100 -b $seed -N 1




