#' Calculate the number of reads aligning to each transcript 
#' author: Arvind Rasi Subramaniam
#' date: Oct 7, 2018

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(org.Sc.sgd.db)
library(GenomicAlignments)
library(glue)
library(tidyverse)
library(plyranges)

# extract genome and annotations
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
orgdb <- org.Sc.sgd.db

# get all transcripts
tx <- transcripts(txdb)

# read in all alignments
aln <- readGAlignments(glue("processeddata/{sample}/genome.bam"))

# find overlap of each alignment with transcripts
overlaps <- findOverlaps(aln, tx)

# assign transcript hit as a column
mcols(aln)[queryHits(overlaps), "txHit"] <- subjectHits(overlaps)

#  find number of reads of aligning to each transcript
tx_read_counts <- aln %>% 
  as_tibble() %>% 
  group_by(txHit, strand) %>% 
  count() %>% 
  rename(alncount = n) %>% 
  arrange(desc(alncount)) %>% 
  write_tsv(glue("processeddata/{sample}/tx_read_counts.tsv"))