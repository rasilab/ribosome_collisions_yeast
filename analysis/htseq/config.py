"""Configuration specific to this analysis

    :Author: Arvind Rasi Subramaniam
    :Date: Nov 11, 2018
"""

# import libraries
import os
import regex as re

# get all raw fastq files
fastq_files =  os.listdir('../../data/htseq/')

# extract barcodes from file names
BARCODES = list()
for file in fastq_files:
    match = re.search('_(?P<barcode>[ACTG]{6})_L00.+.fastq.gz', file)
    if match:
        BARCODES.append(match.group('barcode'))

# convert to list of unique barcodes
BARCODES = list(set(BARCODES))
