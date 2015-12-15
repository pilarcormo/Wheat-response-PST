
##Fasta cDNA dumps and Fasta CDS dumps

These files hold the cDNA sequences corresponding to Ensembl Genomes genes.

####FILE NAMES

The files are consistently named following this pattern:
<species>.<assembly>.<eg_version>.<sequence type>.<status>.fa.gz

<species>: The systematic name of the species. 
<assembly>: The assembly build name.
<eg_version>: The version of Ensembl Genomes from which the data was exported.
<sequence type>: cdna for cDNA sequences
<status>

 -  **'cdna.all'** - the super-set of all transcripts resulting from 
     Ensembl known, novel and pseudo gene predictions (see more below).
     
 - **'cdna.abinitio'** - transcripts resulting from 'ab initio' gene prediction 
     algorithms such as SNAP and GENSCAN. In general all 'ab initio' 
     predictions are solely based on the genomic sequence and do not 
     use other experimental evidence. Therefore, not all GENSCAN or SNAP 
     cDNA predictions represent biologically real cDNAs. 
     Consequently, these predictions should be used with care.
     
  -  **'cds.all'** - the super-set of all transcripts resulting from 
     Ensembl known, novel and pseudo gene predictions 


EXAMPLES  (Note: Most species do not sequences for each different <status>)
  for Human:
    Homo_sapiens.NCBI36.cdna.all.fa.gz
      cDNA sequences for all transcripts: known, novel and pseudo
    Homo_sapiens.NCBI36.cdna.abinitio.fa.gz
      cDNA sequences for 'ab-initio' prediction transcripts.

#####Difference between known and novel transcripts

Transcript or protein models that can be mapped to species-specific entries 
in Swiss-Prot, RefSeq or SPTrEMBL are referred to as known genes in Ensembl.  
Those that cannot be mapped are called novel genes (e.g. genes predicted on 
the basis of evidence from closely related species).









