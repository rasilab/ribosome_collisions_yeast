"""Count reads for each Kozak region and gene by deep sequencing
"""

# Import libraries and define analysis specific variables

import sys
import os
import re
import gzip
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# sequence preceding the -4 to -1 of ATG, should be there in all reads
r1_adapter5 = 'GTTTCGAATAAACACACATCTAGAAA'
# -4 to -1 + start codon
kozakstart = '(?P<kozakstart>[ACTG]{7})'
# 9 nt following start codon
r1_geneseq = '(?P<geneseq>[ACTG]{9})'
# look for all the above 3 patterns in the read continuously
r1_pattern = f'{r1_adapter5}{kozakstart}{r1_geneseq}'

# these are the first 9 nt of each gene that identifies
# the gene uniquely
genestarts = {
    'TCTTTATCT': 'pgk1',   # oHP258
    'ACAGAGCAG': 'his3',  # oHP259
    'GTTTCTGAA': 'optmkate2',  # oHP260
    'GTATCGGAG': 'deoptmkate2',  # oHP261
    'GATTACAAG': 'flag',  # oHP47, this is the generic 3xflag start used often
}


# Get the barcode sequence and corresponding files from command line

# get the barcode from the command line
barcode = sys.argv[1]
# test barcode for interactive testing
# barcode = "GCCTAA"

datadir = '../../../data/htseq/initiation_pgk1_no_insert/'
outputdir = '../../../tables/htseq/initiation_pgk1_no_insert/'
# get all fastq files
fqfiles = os.listdir(datadir)
# get files for this barcode and only forward read
r1files = list(filter(lambda x: barcode in x and "R1" in x, fqfiles))


# Iterate through reads and extract initation region and gene

read_counts = dict()
# keep track of total reads in the file
read_counts['total_reads'] = 0
# reads that have initiation mutation and correct gene
read_counts['init_gene_reads'] = 0

init_gene_counts = dict()

# iterate through barcode fwd read files
for file in r1files:
    # open file for iterating
    r1file = FastqGeneralIterator(gzip.open(f'{datadir}{file}', 'rt'))

    # iterate through records one by one
    for r1 in r1file:
        read_counts['total_reads'] += 1

        # go to next read if 5' adapter is not present
        if r1_adapter5 not in r1[1]:
            continue

#         if read_counts['init_gene_reads'] % 1e5 == 0:
#             print(read_counts)

        # get initiation mutation and the gene sequence
        try:
            r1_match = re.search(r1_pattern, r1[1])
            init_mutation = r1_match.group('kozakstart')
            geneseq = r1_match.group('geneseq')
        except AttributeError or KeyError:
            continue

        try:
            gene = genestarts[geneseq]
        except KeyError:
            continue

        #  increase count for this initiation mutation, gene pair
        key = (init_mutation, gene)
        if key in init_gene_counts:
            init_gene_counts[key] += 1
        else:
            # initialize if this is the first time we are seeing it
            init_gene_counts[key] = 1

        read_counts['init_gene_reads'] += 1


# Write counts for each initiation region and gene to output

data = list()
for k, v in init_gene_counts.items():
    pairdata = dict()
    pairdata['init'] = k[0]
    pairdata['gene'] = k[1]
    pairdata['count'] = v
    data.append(pairdata)

data = pd.DataFrame(data).sort_values(by='count', ascending=False)

data.to_csv(f"{outputdir}{barcode}_init_gene_pairs.tsv",
            sep="\t", index=False)
