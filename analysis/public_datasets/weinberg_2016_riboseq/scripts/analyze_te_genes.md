Analyze TE of yeast genes with or without stalls
================
rasi
17 July, 2019

-   [Load libraries](#load-libraries)
-   [Load genome and annotations](#load-genome-and-annotations)
-   [Load RQC stalls for joining with high TE genes](#load-rqc-stalls-for-joining-with-high-te-genes)
-   [Calculate TE as log2 ratio of RPF to RNA RPKM in Weinberg 2016](#calculate-te-as-log2-ratio-of-rpf-to-rna-rpkm-in-weinberg-2016)
-   [Look at high TE genes with potential RQC stalls](#look-at-high-te-genes-with-potential-rqc-stalls)
-   [Plot TE as a function of stall strength](#plot-te-as-a-function-of-stall-strength)
-   [Test if stall-containing genes have lower or higher TE](#test-if-stall-containing-genes-have-lower-or-higher-te)
-   [Avg TE of stall-containing and other genes](#avg-te-of-stall-containing-and-other-genes)
-   [Session Info](#session-info)

Load libraries
==============

``` r
library(topGO)
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=3)
```

Load genome and annotations
===========================

``` r
annotations <- "/fh/fast/subramaniam_a/db/rasi/genomes/yeast/Saccharomyces_cerevisiae/sgd/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff" %>% 
  rtracklayer::readGFF() %>% 
  as_tibble()
```

Load RQC stalls for joining with high TE genes
==============================================

``` r
rqc_stalls <- read_tsv("../../rqc_stalls_in_yeast_orfs/tables/ngrams_annotated.tsv") %>% 
  select(id, pos, ngram_weight, ngram) %>%
  print()
```

    ## # A tibble: 1,251 x 4
    ##    id        pos ngram_weight ngram     
    ##    <chr>   <int>        <int> <chr>     
    ##  1 YHR131C   321           10 RRRRRRRRRR
    ##  2 YIL159W   768           10 PPPPPPPPPP
    ##  3 YNL271C  1238           10 PPPPPPPPPP
    ##  4 YOR019W   712           10 KKKKKKKKKK
    ##  5 YBL091C    40            9 KKKKNKKKKK
    ##  6 YDL146W   473            9 APPPPPPPPP
    ##  7 YDL173W   273            9 KEKKKKKKKK
    ##  8 YDR174W   227            9 KKKKKKKKDK
    ##  9 YDR334W   231            9 KKKRGRKKKK
    ## 10 YHL019C   174            9 KRKDKKKKRK
    ## # ... with 1,241 more rows

Calculate TE as log2 ratio of RPF to RNA RPKM in Weinberg 2016
==============================================================

``` r
te_data <- list.files("../annotations/", pattern = "RPKMs.txt.gz", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(sample = str_extract(file, "RPF|RiboZero")) %>% 
  mutate(data = map(file, . %>% read_tsv(col_names = F))) %>% 
  select(-sno, -file) %>% 
  unnest() %>% 
  magrittr::set_colnames(c("sampletype", "id", "rpkm")) %>% 
  spread(sampletype, rpkm) %>% 
  dplyr::rename(ribo = RPF, rna = RiboZero) %>% 
  mutate(te = log2(ribo / rna)) %>% 
  filter(!is.na(te) & ribo > 5 & rna > 5) %>%
  left_join(annotations %>% select(gene, ID, Note), by = c("id" = "ID")) %>%
  mutate(Note = as.character(Note)) %>% 
  print()
```

    ## # A tibble: 4,507 x 6
    ##    id        rna    ribo      te gene  Note                               
    ##    <chr>   <dbl>   <dbl>   <dbl> <chr> <chr>                              
    ##  1 YAL00…   37.8   18.1  -1.06   TFC3  Largest of six subunits of the RNA…
    ##  2 YAL00…   22.2    5.73 -1.95   VPS8  Membrane-associated protein that i…
    ##  3 YAL00… 4622.  7065.    0.612  EFB1  Translation elongation factor 1 be…
    ##  4 YAL00…  212.   209.   -0.0202 ERP2  Protein that forms a heterotrimeri…
    ##  5 YAL00…   50.0   95.6   0.935  FUN14 Mitochondrial protein of unknown f…
    ##  6 YAL00…   35.0   18.9  -0.888  SPO7  Putative regulatory subunit of Nem…
    ##  7 YAL01…   16.7   11.0  -0.593  MDM10 Subunit of both the ERMES complex …
    ##  8 YAL01…   17.4   10.5  -0.729  SWC3  Protein of unknown function, compo…
    ##  9 YAL01…  568.  1094.    0.945  CYS3  Cystathionine gamma-lyase, catalyz…
    ## 10 YAL01…   25.5   40.7   0.673  DEP1  Transcriptional modulator involved…
    ## # ... with 4,497 more rows

Look at high TE genes with potential RQC stalls
===============================================

``` r
te_stall_data <- te_data %>% 
  left_join(rqc_stalls, by = "id") %>% 
  arrange(desc(te)) %>% 
  select(te, gene, pos, ngram, ngram_weight, everything()) %>% 
  print()
```

    ## # A tibble: 4,507 x 9
    ##       te gene     pos ngram  ngram_weight id       rna   ribo Note        
    ##    <dbl> <chr>  <int> <chr>         <int> <chr>  <dbl>  <dbl> <chr>       
    ##  1  2.89 TMA23    157 KKKKK…            9 YMR2… 7.46e0   55.3 Nucleolar p…
    ##  2  1.75 MFA1      NA <NA>             NA YDR4… 4.74e2 1596.  Mating pher…
    ##  3  1.69 VPS24     NA <NA>             NA YKL0… 1.96e1   63.3 One of four…
    ##  4  1.65 NOP16    215 KRRLL…            7 YER0… 6.71e1  211.  Constituent…
    ##  5  1.65 TMA10     NA <NA>             NA YLR3… 4.35e1  136.  Protein of …
    ##  6  1.59 TIM17     NA <NA>             NA YJL1… 9.38e1  283.  Essential s…
    ##  7  1.56 <NA>      NA <NA>             NA YDL0… 7.31e1  216.  Tail-anchor…
    ##  8  1.52 RPS29A    NA <NA>             NA YLR3… 2.75e3 7897.  Protein com…
    ##  9  1.47 MFA2      NA <NA>             NA YNL1… 6.44e2 1781.  Mating pher…
    ## 10  1.47 CIS3      NA <NA>             NA YJL1… 6.93e2 1918.  Mannose-con…
    ## # ... with 4,497 more rows

Plot TE as a function of stall strength
=======================================

``` r
plot_data <- te_stall_data %>% 
  mutate(ngram_weight = as.factor(if_else(is.na(ngram_weight), 0, 1))) %>% 
  group_by(ngram_weight) %>% 
  mutate(`n` = paste0("N = ", dplyr::n())) %>% 
  ungroup() %>% 
  mutate(ngram_weight = fct_recode(ngram_weight, `No stall` = "0", `Stall` = "1"))
  
plot_data %>% 
  ggplot(aes(x = ngram_weight, y = te, fill = ngram_weight)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  labs(x = "S. cerevisiae genes", y = "Translation efficiency (log2, a.u.)") +
  geom_text(aes(x = ngram_weight, label = n),
            data = plot_data %>% group_by(ngram_weight) %>% slice(1),
            y = -7, size = 2.8) +
  scale_y_continuous(limits = c(-7.2, 3)) +
  scale_fill_manual(values = cbPalette, guide = "none") +
  NULL
```

![](analyze_te_genes_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
ggsave("../figures/distribution_of_translation_efficiency_for_rqc_stall_containing_saccer_genes.pdf")
```

Test if stall-containing genes have lower or higher TE
======================================================

``` r
wilcox.test(te ~ ngram_weight, data = plot_data, alternative = "two.sided")
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  te by ngram_weight
    ## W = 2111700, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

Avg TE of stall-containing and other genes
==========================================

``` r
plot_data %>% 
  group_by(ngram_weight) %>% 
  summarize(te = mean(te), n = dplyr::n()) %>% 
  knitr::kable()
```

| ngram\_weight |          te|     n|
|:--------------|-----------:|-----:|
| No stall      |  -0.3292908|  3477|
| Stall         |  -0.5665720|  1030|

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
    ##  [1] bindrcpp_0.2.2          rasilabRtemplates_0.1.0
    ##  [3] forcats_0.3.0           stringr_1.3.1          
    ##  [5] dplyr_0.7.6             purrr_0.2.5            
    ##  [7] readr_1.1.1             tidyr_0.8.1            
    ##  [9] tibble_1.4.2            ggplot2_3.0.0          
    ## [11] tidyverse_1.2.1         topGO_2.32.0           
    ## [13] SparseM_1.77            GO.db_3.6.0            
    ## [15] AnnotationDbi_1.42.1    IRanges_2.14.10        
    ## [17] S4Vectors_0.18.3        Biobase_2.40.0         
    ## [19] graph_1.58.0            BiocGenerics_0.26.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137                bitops_1.0-6               
    ##  [3] matrixStats_0.54.0          lubridate_1.7.4            
    ##  [5] bit64_0.9-7                 httr_1.3.1                 
    ##  [7] rprojroot_1.3-2             GenomeInfoDb_1.16.0        
    ##  [9] tools_3.5.1                 backports_1.1.2            
    ## [11] utf8_1.1.4                  R6_2.2.2                   
    ## [13] DBI_1.0.0                   lazyeval_0.2.1             
    ## [15] colorspace_1.3-2            withr_2.1.2                
    ## [17] tidyselect_0.2.4            bit_1.1-14                 
    ## [19] compiler_3.5.1              cli_1.0.0                  
    ## [21] rvest_0.3.2                 xml2_1.2.0                 
    ## [23] DelayedArray_0.6.2          labeling_0.3               
    ## [25] rtracklayer_1.40.3          scales_0.5.0               
    ## [27] digest_0.6.15               Rsamtools_1.32.2           
    ## [29] rmarkdown_1.10              XVector_0.20.0             
    ## [31] pkgconfig_2.0.1             htmltools_0.3.6            
    ## [33] highr_0.7                   rlang_0.2.1                
    ## [35] readxl_1.1.0                rstudioapi_0.7             
    ## [37] RSQLite_2.1.1               bindr_0.1.1                
    ## [39] jsonlite_1.5                BiocParallel_1.14.2        
    ## [41] RCurl_1.95-4.11             magrittr_1.5               
    ## [43] GenomeInfoDbData_1.1.0      Matrix_1.2-14              
    ## [45] fansi_0.2.3                 Rcpp_0.12.18               
    ## [47] munsell_0.5.0               stringi_1.2.4              
    ## [49] yaml_2.2.0                  SummarizedExperiment_1.10.1
    ## [51] zlibbioc_1.26.0             plyr_1.8.4                 
    ## [53] grid_3.5.1                  blob_1.1.1                 
    ## [55] crayon_1.3.4                lattice_0.20-35            
    ## [57] Biostrings_2.48.0           haven_1.1.2                
    ## [59] hms_0.4.2                   knitr_1.20                 
    ## [61] pillar_1.3.0                GenomicRanges_1.32.6       
    ## [63] XML_3.98-1.12               glue_1.3.0                 
    ## [65] evaluate_0.11               modelr_0.1.2               
    ## [67] cellranger_1.1.0            gtable_0.2.0               
    ## [69] assertthat_0.2.0            broom_0.5.0                
    ## [71] GenomicAlignments_1.16.0    memoise_1.1.0
