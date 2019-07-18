Plot ribosome density around RQC stall in endogenous genes
================
rasi
17 July, 2019

-   [Load libraries](#load-libraries)
-   [Load genome and annotations](#load-genome-and-annotations)
-   [Load RQC stalls for joining with high TE genes](#load-rqc-stalls-for-joining-with-high-te-genes)
-   [Convert RQC stalls to genomic coordinates](#convert-rqc-stalls-to-genomic-coordinates)
-   [Load the alignments](#load-the-alignments)
-   [Trim the alignments to the P-site and calculate coverage separately for + and - strands](#trim-the-alignments-to-the-p-site-and-calculate-coverage-separately-for-and---strands)
-   [Load pre-computed coverage](#load-pre-computed-coverage)
-   [Expand each stall to 61 nt width](#expand-each-stall-to-61-nt-width)
-   [Get the coverage of the 61 nt window arounde each stall](#get-the-coverage-of-the-61-nt-window-arounde-each-stall)
-   [Normalize the cvg within each stall and threshold to stalls with mean cvg &gt;= 1](#normalize-the-cvg-within-each-stall-and-threshold-to-stalls-with-mean-cvg-1)
-   [Plot the mean and standard deviation of normalized cvg around stalls](#plot-the-mean-and-standard-deviation-of-normalized-cvg-around-stalls)
-   [Session Info](#session-info)

Load libraries
==============

``` r
library(GenomicAlignments)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)
library(plyranges)
library(rasilabRtemplates)
```

Load genome and annotations
===========================

``` r
genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
annotations <- "/fh/fast/subramaniam_a/db/rasi/genomes/yeast/Saccharomyces_cerevisiae/sgd/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff" %>% 
  rtracklayer::readGFF() %>% 
  as_tibble()

tx <- annotations %>% 
  GRanges() %>% 
  filter(type == "CDS") %>% 
  select(Name) %>% 
  split(.$Name) %>% 
  print()
```

    ## GRangesList object of length 6717:
    ## $Q0010 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames    ranges strand |        Name
    ##          <Rle> <IRanges>  <Rle> | <character>
    ##   [1]  chrMito 3952-4338      + |       Q0010
    ## 
    ## $Q0017 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames    ranges strand |  Name
    ##   [1]  chrMito 4254-4415      + | Q0017
    ## 
    ## $Q0032 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames      ranges strand |  Name
    ##   [1]  chrMito 11667-11957      + | Q0032
    ## 
    ## ...
    ## <6714 more elements>
    ## -------
    ## seqinfo: 18 sequences from an unspecified genome; no seqlengths

Load RQC stalls for joining with high TE genes
==============================================

``` r
rqc_stalls <- read_tsv("../../rqc_stalls_in_yeast_orfs/tables/ngrams_annotated.tsv") %>% 
  bind_rows(read_tsv("../../rqc_stalls_in_yeast_orfs/tables/ngram_control_annotated.tsv")) %>%
  mutate(seqname = id, start = pos*3 + 1) %>% 
  mutate(end = start) %>%
  select(seqname, start, end, stall, id) %>% 
  GRanges() %>% 
  mutate(stall = if_else(stall %in% c("KR", "P"), "KPR", stall)) %>% 
  print()
```

    ## GRanges object with 2803 ranges and 2 metadata columns:
    ##          seqnames    ranges strand |       stall          id
    ##             <Rle> <IRanges>  <Rle> | <character> <character>
    ##      [1]  YHR131C       964      * |         KPR     YHR131C
    ##      [2]  YIL159W      2305      * |         KPR     YIL159W
    ##      [3]  YNL271C      3715      * |         KPR     YNL271C
    ##      [4]  YOR019W      2137      * |         KPR     YOR019W
    ##      [5]  YBL091C       121      * |         KPR     YBL091C
    ##      ...      ...       ...    ... .         ...         ...
    ##   [2799]  YPR164W       436      * |          DE     YPR164W
    ##   [2800]  YPR169W      1321      * |          DE     YPR169W
    ##   [2801]  YPR173C       304      * |          DE     YPR173C
    ##   [2802]  YPR180W       601      * |          DE     YPR180W
    ##   [2803]  YPR204W      2281      * |          DE     YPR204W
    ##   -------
    ##   seqinfo: 2244 sequences from an unspecified genome; no seqlengths

Convert RQC stalls to genomic coordinates
=========================================

``` r
rqc_stalls_coords <- mapFromTranscripts(rqc_stalls, tx) %>% 
  # get rid of mitochondrial sequence
  filter(seqnames != "chrMito") %>% 
  mutate(id = seqnames(rqc_stalls)[xHits], stall = rqc_stalls$stall[xHits]) %>% 
  select(-xHits, -transcriptsHits)

# check that the mapping was done correctly
rqc_stalls_coords %>% 
  anchor_5p() %>% 
  stretch(29) %>% 
  getSeq(genome, .) %>% 
  translate()
```

    ##   A AAStringSet instance of length 2802
    ##        width seq
    ##    [1]    10 RRRRRRRRRR
    ##    [2]    10 PPPPPPPPPP
    ##    [3]    10 PPPPPPPPPP
    ##    [4]    10 KKKKKKKKKK
    ##    [5]    10 KKKKNKKKKK
    ##    ...   ... ...
    ## [2798]    10 DDSDIEALDD
    ## [2799]    10 EELLDELDKD
    ## [2800]    10 EEGEDNGGED
    ## [2801]    10 DEEDEKKTYE
    ## [2802]    10 DEEYKEYLED

Load the alignments
===================

We do not run the codecell below after the first time to save time.

``` r
aln <- readGAlignments("../processeddata/mono/accepted_hits.bam") %>% 
  print()
```

Trim the alignments to the P-site and calculate coverage separately for + and - strands
=======================================================================================

We do not run the codecell below after the first time to save time.

``` r
cvg_plus <- aln[strand(aln) == "+"] %>% 
  qnarrow(start = 13, width = 1) %>% 
  coverage() %>% 
  print()

cvg_minus <- aln[strand(aln) == "-"] %>% 
  qnarrow(start = qwidth(.) - 12, width = 1) %>% 
  coverage() %>% 
  print()

rtracklayer::export.bw(cvg_plus, "../processeddata/mono/cvg_plus.bw")
rtracklayer::export.bw(cvg_minus, "../processeddata/mono/cvg_minus.bw")
```

Load pre-computed coverage
==========================

``` r
cvg_plus <- rtracklayer::import.bw("../processeddata/mono/cvg_plus.bw") %>% 
  coverage(weight = "score")
cvg_minus <- rtracklayer::import.bw("../processeddata/mono/cvg_minus.bw") %>% 
  coverage(weight = "score")
```

Expand each stall to 61 nt width
================================

``` r
rqc_flank <- rqc_stalls_coords %>% 
  anchor_5p() %>% 
  stretch(300) %>% 
  shift_upstream(150) %>% 
  mutate(uid = paste0(id, "_", stall)) %>% 
  print()
```

    ## GRanges object with 2802 ranges and 3 metadata columns:
    ##          seqnames        ranges strand |      id       stall         uid
    ##             <Rle>     <IRanges>  <Rle> |   <Rle> <character> <character>
    ##      [1]  chrVIII 366779-367079      - | YHR131C         KPR YHR131C_KPR
    ##      [2]    chrIX   43979-44279      + | YIL159W         KPR YIL159W_KPR
    ##      [3]   chrXIV 131519-131819      - | YNL271C         KPR YNL271C_KPR
    ##      [4]    chrXV 370113-370413      + | YOR019W         KPR YOR019W_KPR
    ##      [5]    chrII   48358-48658      - | YBL091C         KPR YBL091C_KPR
    ##      ...      ...           ...    ... .     ...         ...         ...
    ##   [2798]   chrXVI 870988-871288      + | YPR164W          DE  YPR164W_DE
    ##   [2799]   chrXVI 879860-880160      + | YPR169W          DE  YPR169W_DE
    ##   [2800]   chrXVI 887384-887684      - | YPR173C          DE  YPR173C_DE
    ##   [2801]   chrXVI 896411-896711      + | YPR180W          DE  YPR180W_DE
    ##   [2802]   chrXVI 946733-947033      + | YPR204W          DE  YPR204W_DE
    ##   -------
    ##   seqinfo: 18 sequences from an unspecified genome; no seqlengths

``` r
rqc_flank_plus <- filter(rqc_flank, strand == "+")
rqc_flank_minus <- filter(rqc_flank, strand == "-")
```

Get the coverage of the 61 nt window arounde each stall
=======================================================

``` r
stall_cvg_plus <- cvg_plus[rqc_flank_plus] %>% 
  setNames(rqc_flank_plus$uid) %>% 
  GRanges() %>% 
  as_tibble()

stall_cvg_minus <- cvg_plus[rqc_flank_minus] %>% 
  setNames(rqc_flank_minus$uid) %>% 
  GRanges() %>% 
  as_tibble()

stall_cvg <- bind_rows(stall_cvg_plus, stall_cvg_minus) %>% 
  # create a sequence from start to stop for each range
  mutate(pos = map2(start, end, function(x, y) seq(from = x, to = y))) %>% 
  # expand each range to equal its length
  unnest()  %>% 
  # mutate and unnest to create a single pos for each location
  mutate(start = pos, end = pos) %>% 
  select(-pos, -width) %>% 
  print()
```

    ## # A tibble: 843,402 x 5
    ##    seqnames    start   end strand score
    ##    <chr>       <int> <int> <fct>  <dbl>
    ##  1 YIL159W_KPR     1     1 *          0
    ##  2 YIL159W_KPR     2     2 *          0
    ##  3 YIL159W_KPR     3     3 *          0
    ##  4 YIL159W_KPR     4     4 *          0
    ##  5 YIL159W_KPR     5     5 *          0
    ##  6 YIL159W_KPR     6     6 *          0
    ##  7 YIL159W_KPR     7     7 *          1
    ##  8 YIL159W_KPR     8     8 *          0
    ##  9 YIL159W_KPR     9     9 *          1
    ## 10 YIL159W_KPR    10    10 *          0
    ## # ... with 843,392 more rows

Normalize the cvg within each stall and threshold to stalls with mean cvg &gt;= 1
=================================================================================

``` r
norm_cvg <- stall_cvg %>% 
  group_by(seqnames) %>% 
  mutate(mean_score = mean(score)) %>% 
  mutate(norm_score = score / mean_score) %>%
  ungroup() %>% 
  filter(mean_score >= 1) %>% 
  print()
```

    ## # A tibble: 280,532 x 7
    ##    seqnames    start   end strand score mean_score norm_score
    ##    <chr>       <int> <int> <fct>  <dbl>      <dbl>      <dbl>
    ##  1 YDL173W_KPR     1     1 *         13       4.92      2.64 
    ##  2 YDL173W_KPR     2     2 *          3       4.92      0.610
    ##  3 YDL173W_KPR     3     3 *          8       4.92      1.63 
    ##  4 YDL173W_KPR     4     4 *         37       4.92      7.52 
    ##  5 YDL173W_KPR     5     5 *          5       4.92      1.02 
    ##  6 YDL173W_KPR     6     6 *         10       4.92      2.03 
    ##  7 YDL173W_KPR     7     7 *         18       4.92      3.66 
    ##  8 YDL173W_KPR     8     8 *          1       4.92      0.203
    ##  9 YDL173W_KPR     9     9 *          3       4.92      0.610
    ## 10 YDL173W_KPR    10    10 *          9       4.92      1.83 
    ## # ... with 280,522 more rows

Plot the mean and standard deviation of normalized cvg around stalls
====================================================================

``` r
plot_data <- norm_cvg %>% 
  mutate(stall = str_extract(seqnames, "[^_]+$")) %>% 
  group_by(start, stall) %>% 
  summarize(mean_density = mean(norm_score), sd_density = sd(norm_score))

plot_data %>% 
  ggplot(aes(x = (start - 150), y = mean_density, 
             ymin = mean_density - sd_density, ymax = mean_density + sd_density)) +
  facet_wrap(~ stall, ncol = 3, scales = "free") +
  scale_y_continuous(limits = c(0, 4)) +
  scale_x_continuous(limits = c(-60, 60)) +
  labs(x = "Distance from endogenous gene stall (nt)", y = "Mean ribosome density (a.u.)") +
  geom_line(size = 0.5) 
```

![](plot_ribo_density_around_rqc_stalls_files/figure-markdown_github/plot_stall_cvg-1.png)

``` r
ggsave("../figures/ribosome_density_around_rqc_stalls_and_controls.pdf")
```

Session Info
============

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /home/rasi/lib/R-3.5.1/lib/R/lib/libRblas.so
    ## LAPACK: /home/rasi/lib/R-3.5.1/lib/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices datasets  utils    
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2              rasilabRtemplates_0.1.0    
    ##  [3] plyranges_1.0.3             forcats_0.3.0              
    ##  [5] stringr_1.3.1               dplyr_0.7.6                
    ##  [7] purrr_0.2.5                 readr_1.1.1                
    ##  [9] tidyr_0.8.1                 tibble_1.4.2               
    ## [11] ggplot2_3.0.0               tidyverse_1.2.1            
    ## [13] GenomicFeatures_1.32.0      AnnotationDbi_1.42.1       
    ## [15] GenomicAlignments_1.16.0    Rsamtools_1.32.2           
    ## [17] Biostrings_2.48.0           XVector_0.20.0             
    ## [19] SummarizedExperiment_1.10.1 DelayedArray_0.6.2         
    ## [21] BiocParallel_1.14.2         matrixStats_0.54.0         
    ## [23] Biobase_2.40.0              GenomicRanges_1.32.6       
    ## [25] GenomeInfoDb_1.16.0         IRanges_2.14.10            
    ## [27] S4Vectors_0.18.3            BiocGenerics_0.26.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137                           
    ##  [2] bitops_1.0-6                           
    ##  [3] lubridate_1.7.4                        
    ##  [4] bit64_0.9-7                            
    ##  [5] progress_1.2.0                         
    ##  [6] httr_1.3.1                             
    ##  [7] rprojroot_1.3-2                        
    ##  [8] BSgenome.Scerevisiae.UCSC.sacCer3_1.4.0
    ##  [9] tools_3.5.1                            
    ## [10] backports_1.1.2                        
    ## [11] utf8_1.1.4                             
    ## [12] R6_2.2.2                               
    ## [13] DBI_1.0.0                              
    ## [14] lazyeval_0.2.1                         
    ## [15] colorspace_1.3-2                       
    ## [16] withr_2.1.2                            
    ## [17] tidyselect_0.2.4                       
    ## [18] prettyunits_1.0.2                      
    ## [19] bit_1.1-14                             
    ## [20] compiler_3.5.1                         
    ## [21] cli_1.0.0                              
    ## [22] rvest_0.3.2                            
    ## [23] xml2_1.2.0                             
    ## [24] labeling_0.3                           
    ## [25] rtracklayer_1.40.3                     
    ## [26] scales_0.5.0                           
    ## [27] digest_0.6.15                          
    ## [28] rmarkdown_1.10                         
    ## [29] pkgconfig_2.0.1                        
    ## [30] htmltools_0.3.6                        
    ## [31] BSgenome_1.48.0                        
    ## [32] rlang_0.2.1                            
    ## [33] readxl_1.1.0                           
    ## [34] rstudioapi_0.7                         
    ## [35] RSQLite_2.1.1                          
    ## [36] bindr_0.1.1                            
    ## [37] jsonlite_1.5                           
    ## [38] RCurl_1.95-4.11                        
    ## [39] magrittr_1.5                           
    ## [40] GenomeInfoDbData_1.1.0                 
    ## [41] Matrix_1.2-14                          
    ## [42] Rcpp_0.12.18                           
    ## [43] munsell_0.5.0                          
    ## [44] fansi_0.2.3                            
    ## [45] stringi_1.2.4                          
    ## [46] yaml_2.2.0                             
    ## [47] zlibbioc_1.26.0                        
    ## [48] plyr_1.8.4                             
    ## [49] grid_3.5.1                             
    ## [50] blob_1.1.1                             
    ## [51] crayon_1.3.4                           
    ## [52] lattice_0.20-35                        
    ## [53] haven_1.1.2                            
    ## [54] hms_0.4.2                              
    ## [55] knitr_1.20                             
    ## [56] pillar_1.3.0                           
    ## [57] biomaRt_2.36.1                         
    ## [58] XML_3.98-1.12                          
    ## [59] glue_1.3.0                             
    ## [60] evaluate_0.11                          
    ## [61] modelr_0.1.2                           
    ## [62] cellranger_1.1.0                       
    ## [63] gtable_0.2.0                           
    ## [64] assertthat_0.2.0                       
    ## [65] broom_0.5.0                            
    ## [66] memoise_1.1.0
