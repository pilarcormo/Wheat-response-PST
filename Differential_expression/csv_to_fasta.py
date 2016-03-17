
##download DE expressed genes with biopython interface to entrez
##usage: python listeria_entrez.py someone@somewhere.someaddress -replaceoutliers-results.csv out.fasta


from Bio import Entrez, SeqIO
import sys

# command line parameters
email = sys.argv[1]
csv = sys.argv[2]
out = sys.argv[3]

# initialize variables
Entrez.email = email
gene_ids = []
database = 'protein'
gene_sequences = []

# extract gene ids from table
for line in open(csv):
    csv_row = line.replace('\"', '').split(',')[0]
    if csv_row.strip() != '':
        gene_ids.append(csv_row)

for gene in gene_ids:
    # use the esearch utility to query the ncbi protein database with gene id
    # return the first sequence
    handle = Entrez.esearch(db = database, retmax = 1, term = gene)
    record = Entrez.read(handle)
    handle.close()
    # try to read search results, throw exception if no results
    try:
        # use the efetch utility to get a record of the esearch results
        fetch_handle = Entrez.efetch(db = database, id = record['IdList'][0], rettype = 'fasta', retmode='text')
        fetch_record = SeqIO.read(fetch_handle, 'fasta')
        fetch_handle.close()
        # add sequence to list
        gene_sequences.append(fetch_record)
        print('Gene: ' + gene + ' was successfully downloaded from the NCBI protein database')
    except IndexError:
        print('Gene: ' + gene + ' was not found in the NCBI protein database')

# open out file and write sequences to file
output_handle = open(out, 'w')
SeqIO.write(gene_sequences, output_handle, 'fasta')
output_handle.close()