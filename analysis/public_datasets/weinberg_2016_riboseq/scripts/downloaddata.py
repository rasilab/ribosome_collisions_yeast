#!/fh/fast/subramaniam_a/user/rasi/virtualenv/default2/bin/python
# -*- coding: utf-8 -*-
'''Download reads in .sra format and convert it to .fastq format
'''

# python modules------------------------------------------------------
from __future__ import print_function
import sys
from datetime import datetime
import os
import subprocess as sp
import pandas as pd

# set up sample specific variables------------------------------------
sample = sys.argv[1]
# read parameters for this sample
sampleinfo = pd.read_table(
    'annotations/srafiles.txt', index_col=0).ix[sample]
# create an output folder for this sample
outputdir = 'processeddata/' + sample + '/'

# quit this script if fastq file is already present in rawdata folder
fastqfile = sample + '.fastq.gz'
if fastqfile in os.listdir('rawdata'):
    quit()

# filename to write processing stats
logfile = open(outputdir + 'processingstats.txt', 'a')

# download sra file
sp.check_output([
    'wget -q -r -nH -nd -np -R index.html* ' +
    sampleinfo['srapath'] +
    ' -O rawdata/' + sample + '.sra'
], shell=True)
logfile.write("Downloaded sra file {0}\n".format(
    str(datetime.now())))

# download fastqfile
sp.check_output([
    'fastq-dump',
    '--gzip',
    '-O',
    'rawdata/',
    'rawdata/' + sample + '.sra',
])
logfile.write("Converted sra to fastq {0}\n".format(
    str(datetime.now())))
logfile.close()
