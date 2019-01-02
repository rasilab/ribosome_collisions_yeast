#!/fh/fast/subramaniam_a/user/rasi/virtualenv/py362/bin/python

"""Calculate the number of reads aligning to each barcode
@author Arvind Subramaniam
@date 11 Nov 2018
"""

## Import libraries and define variables
import sys
import os
import regex as re
from Bio import Seq
import gzip
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# We are PCRing from the cyc1 terminator towards YFP
# This is the sequence preceding barcode and is part of the terminator
r1_adapter5 = 'GTGACATAACTAATTACATGA'
# This is the barcode we are looking for
r1_barcode = '(?P<r1_barcode>[ACTG]{8})'
# This is the sequence following barcode and is the restriction site + stop
# codon
r1_adapter3 = 'GTCGACTTA'
r1_pattern = re.compile(f'{r1_adapter5}{r1_barcode}{r1_adapter3}')

# folder with data and annotations
data_folder = '../../data/htseq'

# Get barcode and raw files

# get the barcode from the command line
barcode = sys.argv[1]
# get all fastq files
fqfiles = os.listdir(data_folder)
# get files for this barcode and only forward read
r1files = list(filter(lambda x: barcode in x  and 'R1' in x, fqfiles))

# we count only the barcodes that were part of our original cloning
r1_barcodes = pd.read_table(f'{data_folder}/barcode_annotations.tsv')['barcode'].tolist()

## Get read counts for each initiation-codon-r1_barcode combination
read_counts = dict()
# keep track of total paired reads in the file
read_counts['total'] = 0
# reads that have the 5' adapter in both fwd and rev reads
read_counts['adapter'] = 0
# reads that have initiation mutation, codon mutation, r1_barcode
read_counts['useful'] = 0

# This dict is used to store the number of reads for each combination
combo_counts = dict()

# iterate through barcode fwd read files
for file in r1files:
    # open file for iterating
    r1file = FastqGeneralIterator(gzip.open(f'{data_folder}/{file}', 'rt'))
    
    # iterate through records one by one
    for r1 in r1file:
        
        read_counts['total'] += 1
        
        # go to next read if 5' adapter is not present in Read 1, or if
        # go to next read if 5' adapter is not present in Read 1
        if r1_adapter5 not in r1[1]:
            continue

        read_counts['adapter'] += 1
        if read_counts['adapter'] % 1e5 == 0:
            print(read_counts)
           
        # get barcode from the read
        try:
            r1_match = r1_pattern.search(r1[1])
            r1_barcode = r1_match.group('r1_barcode')
        except AttributeError:
            continue

        # include only barcodes that we cloned
        if r1_barcode not in r1_barcodes:
            continue
    
        #  increase count for this r1_barcode
        key = r1_barcode
        if key in combo_counts:
            combo_counts[key] += 1
        else:
            # initialize if this is the first time we are seeing it
            combo_counts[key] = 1
            
        read_counts['useful'] += 1

## Write combination counts counts to tsv
# first convert it to a Dataframe through a list
data = list()
for k,v in combo_counts.items():
    combodata = dict()
    combodata['r1_barcode'] = k
    combodata['count'] = v
    data.append(combodata)
(
    pd.DataFrame(data).
    sort_values(by = ['r1_barcode', 'count'], ascending=False).
    reset_index().
    to_csv(f'tables/{barcode}_r1barcode_counts.tsv',
           sep='\t', 
           index=False)
)

## Write alignment statistics to tsv
alignment_stats = pd.DataFrame(pd.Series(read_counts).
                               sort_values(ascending=False)).reset_index()
alignment_stats.columns = ['type', 'count']
alignment_stats.to_csv(f'tables/{barcode}_read_counts.tsv',
                       sep='\t',
                       index=False)
