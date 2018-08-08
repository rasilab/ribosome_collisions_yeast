"""Calculate number of reads aligning to each Kozak and stall sequence
"""

# Import libraries and define variables
import sys
import os
import re
from Bio import Seq
import gzip
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

datadir = '../../../data/htseq/initiation_pgk1_aag_ccg/'
outputdir = '../../../tables/htseq/initiation_pgk1_aag_ccg/'

r2_adapter5 = 'GTTTCGAATAAACACACATCTAGAAA'
kozakstart = '(?P<kozakstart>[ACTG]{7})'
r2_pattern = f'{r2_adapter5}{kozakstart}'

r1_adapter5 = 'CAACACCAGTGAACAATTCTTCACCCTTGGATCC'
codonmutation = '(?P<codonmutation>[ACTG]{12})'
r1_pattern = f'{r1_adapter5}{codonmutation}'


# Get barcode and raw files

# get the barcode from the command line
barcode = sys.argv[1]
# test barcode for interactive testing
# barcode = "CGTGAT"

# get all fastq files
fqfiles = os.listdir(datadir)
# get files for this barcode and only forward read
r1files = list(filter(lambda x: barcode in x and "R1" in x, fqfiles))


# Get read counts for each initiation-codon pair
read_counts = dict()
# keep track of total paired reads in the file
read_counts['total_reads'] = 0
# reads that have the 5' adapter in both fwd and rev reads
read_counts['both_adapter_reads'] = 0
# reads that have initiation mutation and codon mutation
read_counts['init_codon_reads'] = 0

init_codon_pair_counts = dict()

# iterate through barcode fwd read files
for file in r1files:
    # open file for iterating
    r1file = FastqGeneralIterator(gzip.open(f'{datadir}{file}', 'rt'))
    # open corresponding reverse read file
    r2filename = file.replace("R1", "R2")
    r2file = FastqGeneralIterator(gzip.open(f'{datadir}{r2filename}', 'rt'))

    # iterate through records one by one
    for r1, r2 in zip(r1file, r2file):
        # check that the ids are the same
        assert(r1[0].split()[0] == r2[0].split()[0])

        read_counts['total_reads'] += 1

        # go to next read if both 5' adapters are not present
        if r1_adapter5 not in r1[1] or r2_adapter5 not in r2[1]:
            continue
        read_counts['both_adapter_reads'] += 1
#         if read_counts['both_adapter_reads'] % 1e5 == 0:
#             print(read_counts)

        # get initiation mutation and the reverse complement codon mutation
        try:
            init_mutation = re.search(r2_pattern, r2[1]).group('kozakstart')
            codon_mutation = re.search(r1_pattern, r1[1]).group('codonmutation')
        except AttributeError:
            continue

        #  increase count for this initiation mutation, codon pair
        key = (init_mutation, codon_mutation)
        if key in init_codon_pair_counts:
            init_codon_pair_counts[key] += 1
        else:
            # initialize if this is the first time we are seeing it
            init_codon_pair_counts[key] = 1

        read_counts['init_codon_reads'] += 1


# Write init-codon pair counts to tsv
data = list()
for k, v in init_codon_pair_counts.items():
    pairdata = dict()
    pairdata['init'] = k[0]
    pairdata['codon'] = Seq.Seq(k[1]).reverse_complement().__str__()
    pairdata['count'] = v
    data.append(pairdata)

(
    pd.DataFrame(data).
    sort_values(by="count", ascending=False).
    reset_index().
    to_csv(f"{outputdir}{barcode}_init_codon_pairs.tsv", sep="\t", index=False)
)

# Write alignment statistics to tsv
(
    pd.Series(read_counts).
    to_csv(f"{outputdir}{barcode}_read_counts.tsv", sep="\t", index=False)
)
