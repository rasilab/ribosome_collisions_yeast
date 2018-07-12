"""Identify RQC-triggering patches in S.cerevisiae proteins
"""

from Bio import SeqIO
import re
import pandas as pd

# look for ngram/stretches of K / R  or P of this length
ngram_length = 10
# the min number of K / R residues in ngrams
ngram_weight_cutoff = 6

ngrams = dict()
# These are patches that are used as controls for analysis
ngram_controls = dict()
# iterate through coding sequences
for rec in SeqIO.parse("../annotations/cds.fa", "fasta"):
    # translate to proteins
    translation = str(rec.seq.translate())
    ngrams[rec.id] = list()
    ngram_controls[rec.id] = list()
    # move through protein sequence
    for i in range(len(translation)):
        # look at stretches of ngram length
        ngram = translation[i:i+ngram_length]
        # look at fraction of K / R residues in this ngram
        ngram_weight = (ngram.count("K") + ngram.count("R"))
        # wirte the sequence, weight, and location if it passes the cutoff
        # above
        if ngram_weight >= ngram_weight_cutoff:
            ngrams[rec.id].append({'ngram': ngram, 'ngram_weight': ngram_weight, "pos": i, "stall": "KR"})
        # look at fraction of P residues in this ngram
        ngram_weight = ngram.count("P") 
        # wirte the sequence, weight, and location if it passes the cutoff
        # above
        if ngram_weight >= ngram_weight_cutoff:
            ngrams[rec.id].append({'ngram': ngram, 'ngram_weight': ngram_weight, "pos": i, "stall": "P"})
        # look at fraction of D/E residues in this ngram
        ngram_weight = (ngram.count("D") + ngram.count("E")) 
        # wirte the sequence, weight, and location if it passes the cutoff
        # above
        if ngram_weight >= ngram_weight_cutoff:
            ngram_controls[rec.id].append({'ngram': ngram, 'ngram_weight': ngram_weight, "pos": i, "stall": "DE"})
    # convert to a dataframe
    ngrams[rec.id] = pd.DataFrame(ngrams[rec.id])
    ngram_controls[rec.id] = pd.DataFrame(ngram_controls[rec.id])

# combine ngrams for all proteins and write it to a tsv file
(
pd.concat(ngrams, names = ["id", "index"])
 .reset_index()
 .drop(["index"], axis=1)
 .to_csv("../tables/ngrams.tsv", index=False, sep="\t")
)
(
pd.concat(ngram_controls, names = ["id", "index"])
 .reset_index()
 .drop(["index"], axis=1)
 .to_csv("../tables/ngram_controls.tsv", index=False, sep="\t")
)
