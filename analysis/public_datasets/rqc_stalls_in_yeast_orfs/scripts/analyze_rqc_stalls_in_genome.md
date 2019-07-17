Identify putative stalls in S. cerevisiae mRNAs that trigger RQC
================
rasi
17 July, 2019

-   [Load libraries and define analysis-specific parameters](#load-libraries-and-define-analysis-specific-parameters)
-   [Get CDS annotations and sequence](#get-cds-annotations-and-sequence)
-   [Load precomputed ngrams](#load-precomputed-ngrams)
-   [Extract annotated RQC sequences of all genes and those for experimental verification](#extract-annotated-rqc-sequences-of-all-genes-and-those-for-experimental-verification)
-   [Extract the RQC controls for plotting ribosome density](#extract-the-rqc-controls-for-plotting-ribosome-density)

This script writes the tables &lt;../tables/ngrams\_annotated.tsv&gt; and &lt;../tables/ngram\_control\_annotated.tsv&gt; that is used by other analysis scripts in the manuscript.

Load libraries and define analysis-specific parameters
======================================================

``` r
library(rasilabRtemplates)
library(biobroom)
library(Biostrings)
library(GenomicFeatures)
library(glue)
library(tidyverse)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
annotations <- "/fh/fast/subramaniam_a/db/rasi/genomes/yeast/Saccharomyces_cerevisiae/sgd/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff" %>% 
  rtracklayer::readGFF()
```

Get CDS annotations and sequence
================================

``` r
cds <- annotations %>% 
  as_tibble() %>% 
  filter(type == "CDS" & !orf_classification == "Dubious") %>% 
  mutate(seqid = if_else(seqid == "chrMito", "chrM", as.character(seqid))) %>% 
  filter(str_detect(seqid, "^chr")) %>% 
  mutate(Note = as.character(Note)) %>% 
  GRanges()

cds %>% 
  split(.$Name) %>% 
  extractTranscriptSeqs(genome, .) %>% 
  writeXStringSet("../annotations/cds.fa")
```

Load precomputed ngrams
=======================

This is done in &lt;count\_rqc\_residues.py&gt;.

``` r
ngrams <- read_tsv("../tables/ngrams.tsv") %>% 
  print()
```

    ## # A tibble: 6,093 x 5
    ##    id      ngram      ngram_weight   pos stall
    ##    <chr>   <chr>             <dbl> <dbl> <chr>
    ##  1 Q0070   KDKNKNKNKK            6   329 KR   
    ##  2 YAL001C KRKRML*RMK            6   292 KR   
    ##  3 YAL011W RTRTRSRRGK            6    20 KR   
    ##  4 YAL011W TRTRSRRGKR            6    21 KR   
    ##  5 YAL011W RTRSRRGKRG            6    22 KR   
    ##  6 YAL011W TRSRRGKRGR            6    23 KR   
    ##  7 YAL011W RSRRGKRGRD            6    24 KR   
    ##  8 YAL011W LKKKEKELKR            6   188 KR   
    ##  9 YAL011W KKKEKELKRK            7   189 KR   
    ## 10 YAL011W KKEKELKRKN            6   190 KR   
    ## # ... with 6,083 more rows

Extract annotated RQC sequences of all genes and those for experimental verification
====================================================================================

``` r
ngrams_annotated <- ngrams %>% 
  left_join(cds %>% tidy() %>% dplyr::select(Name, gene, Note), by = c("id" = "Name")) %>% 
  group_by(id) %>% 
  # pick the severest stall and the one at the 5' most end
  arrange(desc(ngram_weight), pos) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(desc(ngram_weight)) %>% 
  select(-Note, everything(), Note) %>% 
  write_tsv("../tables/ngrams_annotated.tsv") %>% 
  print()
```

    ## # A tibble: 1,251 x 7
    ##    id     ngram   ngram_weight   pos stall gene  Note                     
    ##    <chr>  <chr>          <dbl> <dbl> <chr> <chr> <chr>                    
    ##  1 YHR13… RRRRRR…           10   321 KR    <NA>  Putative protein of unkn…
    ##  2 YIL15… PPPPPP…           10   768 P     BNR1  Formin, nucleates the fo…
    ##  3 YNL27… PPPPPP…           10  1238 P     BNI1  Formin, nucleates the fo…
    ##  4 YOR01… KKKKKK…           10   712 KR    <NA>  Protein of unknown funct…
    ##  5 YBL09… KKKKNK…            9    40 KR    MAP2  Methionine aminopeptidas…
    ##  6 YDL14… APPPPP…            9   473 P     LDB17 Protein involved in the …
    ##  7 YDL17… KEKKKK…            9   273 KR    PAR32 Putative protein of unkn…
    ##  8 YDR17… KKKKKK…            9   227 KR    HMO1  Chromatin associated hig…
    ##  9 YDR33… KKKRGR…            9   231 KR    SWR1  Swi2/Snf2-related ATPase…
    ## 10 YHL01… KRKDKK…            9   174 KR    APM2  Protein of unknown funct…
    ## # ... with 1,241 more rows

Extract the RQC controls for plotting ribosome density
======================================================

``` r
read_tsv("../tables/ngram_controls.tsv") %>% 
  left_join(cds %>% tidy() %>% dplyr::select(Name, gene, Note), by = c("id" = "Name")) %>% 
  group_by(id) %>% 
  # pick the severest stall and the one at the 5' most end
  arrange(desc(ngram_weight), pos) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(desc(ngram_weight)) %>% 
  select(-Note, everything(), Note) %>% 
  write_tsv("../tables/ngram_control_annotated.tsv") %>% 
  print()
```

    ## # A tibble: 1,552 x 7
    ##    id     ngram   ngram_weight   pos stall gene  Note                     
    ##    <chr>  <chr>          <dbl> <dbl> <chr> <chr> <chr>                    
    ##  1 YAL01… DDDDDD…           10    33 DE    SWC3  Protein of unknown funct…
    ##  2 YBL01… EEEEEE…           10   735 DE    SCT1  Glycerol 3-phosphate/dih…
    ##  3 YBR01… EDDDDD…           10   379 DE    KAP1… Transportin or cytosolic…
    ##  4 YBR08… EEEEDD…           10   700 DE    SPT7  Subunit of the SAGA tran…
    ##  5 YBR15… DDDDDD…           10   991 DE    TBS1  Putative protein of unkn…
    ##  6 YBR24… DDEDED…           10    96 DE    ENP1  Protein associated with …
    ##  7 YDL05… DEEDDE…           10  1771 DE    USO1  Essential protein involv…
    ##  8 YDL15… EDEEEE…           10    76 DE    SAS10 Essential subunit of U3-…
    ##  9 YDR02… DDEEEE…           10   741 DE    REG1  Regulatory subunit of ty…
    ## 10 YDR06… DDEEEE…           10   252 DE    DOS2  Protein of unknown funct…
    ## # ... with 1,542 more rows
