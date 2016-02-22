###Variety confirmation

1. Copy the ```Identify_varieties``` directory

2. The reference contigs ```Reference_contigs/contigs_with_markers.fasta```

3. Filter and trim the fastq reads if it's not already done. 
4. Align the reads against the reference_contigs using top_hat
5. Sort and index the BAM files
6. SNP calling
7. Obtain the reference tab file
8. Confirm de variety:

```
ruby Identify_varieties/SNP_markers.rb $m $path
ruby Identify_varieties/decode_scores.rb $m $path
```

Need to specify as arguments: 

- ```$m``` Name of the library 
- ```$path``` path where the library can be found


The pipeline to run it is at
```confirm-variety.sh```



With the output cvs file ```final_scores_<library>.csv``` go to R and run:

```
LIB <- read.csv("final_scores_<library>.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)
varieties <- colSums(LIB[,-1])
df <- data.frame(varieties)
colMax(df)
varieties
```

You will observe something like this:

![image](Screen Shot 2016-02-04 at 15.24.17.png)


Look for the variety that matches the number (335) and that will be the predicted variety -> SOLSTICE