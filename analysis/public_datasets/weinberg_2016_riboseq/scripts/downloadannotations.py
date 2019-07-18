#!/fh/fast/subramaniam_a/user/rasi/virtualenv/default2/bin/python

'''Download annotations from NCBI GEO'''

# libaries------------------------------------------------------------
import subprocess as sp
import gzip
from Bio import Geo
import copy
import pandas as pd
import os
import sys

# samples to be downloaded--------------------------------------------
geoids = {
    # Weinberg et al. Cell Reports 2016
    # Improved Ribosome-Footprint and mRNA Measurements ...
    'GSE53313': {
        'name': 'weinberg_2016',
        'readlength': 40,  # from fastq file
        'names': {  # to make sample names nicer for processing
            'Cerevisiae_RNA': 'totalrz',
            'Cerevisiae_RPF': 'mono',
        }
    }
}

# download annotations------------------------------------------------
samples = dict()
for geoid in geoids:
    url = (
        'ftp://ftp.ncbi.nlm.nih.gov/geo/series/' +
        geoid[:-3] + 'nnn/' +
        geoid + '/soft/' +
        geoid + '_family.soft.gz'
    )
    # download corresponding SOFT data file
    sp.check_output([
        'wget', '-qO', 'annotations/soft.txt.gz', url
    ])
    handle = gzip.open('annotations/soft.txt.gz')
    records = Geo.parse(handle)
    for record in records:
        if record.entity_type == 'SAMPLE':
            # get all supplementary files
            supp_files = [
                record.entity_attributes[x]
                for x in record.entity_attributes
                if x.find('Sample_supplementary_file') != -1]
            # get only the sra file locations
            srafiles = filter(
                lambda x: x.find('sra') != -1, supp_files
            )
            # append sampleset name to samplename
            sample = geoids[geoid]['names'][
                record.entity_attributes['Sample_title']
            ]
            samples[sample] = copy.deepcopy(geoids[geoid])
            samples[sample]['srapath'] = srafiles[0]

# write annotations to console----------------------------------------
samples = pd.DataFrame.from_dict(
    samples, orient='index').drop(['name', 'names'], axis=1)
samples.to_csv('annotations/srafiles.txt', sep='\t')
