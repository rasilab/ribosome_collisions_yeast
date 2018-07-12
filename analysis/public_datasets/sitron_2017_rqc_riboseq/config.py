"""Configuration specific to this analysis

    :Author: Arvind Rasi Subramaniam
    :Datea: Oct 7, 2018
"""

import os
import re

fastq_files =  os.listdir('rawdata/')

RIBO_SAMPLES = list()
MRNA_SAMPLES = list()
for file in fastq_files:
    ribomatch = re.search('([^/]+_mono).fastq.gz', file)
    mrnamatch = re.search('([^/]+_mrna).fastq.gz', file)
    if ribomatch:
        RIBO_SAMPLES.append(ribomatch.groups()[0])
    elif mrnamatch:
        MRNA_SAMPLES.append(mrnamatch.groups()[0])
