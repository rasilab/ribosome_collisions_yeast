Plot ribosome density around RQC stall in endogenous genes
================
rasi
30 July, 2019

-   [Load libraries](#load-libraries)
-   [Load genome and annotations](#load-genome-and-annotations)
-   [Get TE from the the paper](#get-te-from-the-the-paper)
-   [Load RQC stalls for joining with high TE genes](#load-rqc-stalls-for-joining-with-high-te-genes)
-   [Find typical (median) location of stall for use in control](#find-typical-median-location-of-stall-for-use-in-control)
-   [Convert RQC stalls to genomic coordinates](#convert-rqc-stalls-to-genomic-coordinates)
-   [Load the Ribo-seq alignments](#load-the-ribo-seq-alignments)
-   [Trim the Ribo-seq alignments to the P-site and calculate coverage separately for + and - strands](#trim-the-ribo-seq-alignments-to-the-p-site-and-calculate-coverage-separately-for-and---strands)
-   [Load the RNA-seq alignments](#load-the-rna-seq-alignments)
-   [Trim the RNA-seq alignments to the P-site and calculate coverage separately for + and - strands](#trim-the-rna-seq-alignments-to-the-p-site-and-calculate-coverage-separately-for-and---strands)
-   [Load pre-computed coverage](#load-pre-computed-coverage)
-   [Plot calculated TE vs. TE pre-calculated by Weinberg et al.](#plot-calculated-te-vs.-te-pre-calculated-by-weinberg-et-al.)
-   [Look at high TE (paper) genes with potential RQC stalls](#look-at-high-te-paper-genes-with-potential-rqc-stalls)
-   [Plot TE from paper as a function of stall strength](#plot-te-from-paper-as-a-function-of-stall-strength)
-   [Test if stall-containing genes have lower or higher TE from paper](#test-if-stall-containing-genes-have-lower-or-higher-te-from-paper)
-   [Plot difference in TE between stall-containing and remaining genes](#plot-difference-in-te-between-stall-containing-and-remaining-genes)
-   [Look at high TE (calculated) genes with potential RQC stalls](#look-at-high-te-calculated-genes-with-potential-rqc-stalls)
-   [Plot calculated TE as a function of stall strength](#plot-calculated-te-as-a-function-of-stall-strength)
-   [Test if stall-containing genes have lower or higher calculated TE](#test-if-stall-containing-genes-have-lower-or-higher-calculated-te)
-   [Plot difference in TE between stall-containing and remaining genes](#plot-difference-in-te-between-stall-containing-and-remaining-genes-1)
-   [Calculate TE for only the region preceding each stall](#calculate-te-for-only-the-region-preceding-each-stall)
-   [Look at high TE (calculated only preceding stall) genes with potential RQC stalls](#look-at-high-te-calculated-only-preceding-stall-genes-with-potential-rqc-stalls)
-   [Plot calculated TE for region preceding stalls as a function of stall strength](#plot-calculated-te-for-region-preceding-stalls-as-a-function-of-stall-strength)
-   [Test if stall-containing genes have lower or higher calculated TE preceding stall](#test-if-stall-containing-genes-have-lower-or-higher-calculated-te-preceding-stall)
-   [Plot difference in TE between stall-containing and remaining genes](#plot-difference-in-te-between-stall-containing-and-remaining-genes-2)
-   [Session Info](#session-info)
-   [Source data for S5 Fig panel B](#source-data-for-s5-fig-panel-b)

Load libraries
==============

``` r
library(GenomicAlignments)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)
library(plyranges)
library(biobroom)
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
  filter(type == "CDS" & orf_classification == "Verified") %>% 
  select(Name) %>% 
  split(.$Name) %>% 
  print()
```

    ## GRangesList object of length 4896:
    ## $Q0045 
    ## GRanges object with 8 ranges and 1 metadata column:
    ##       seqnames      ranges strand |        Name
    ##          <Rle>   <IRanges>  <Rle> | <character>
    ##   [1]  chrMito 13818-13986      + |       Q0045
    ##   [2]  chrMito 16435-16470      + |       Q0045
    ##   [3]  chrMito 18954-18991      + |       Q0045
    ##   [4]  chrMito 20508-20984      + |       Q0045
    ##   [5]  chrMito 21995-22246      + |       Q0045
    ##   [6]  chrMito 23612-23746      + |       Q0045
    ##   [7]  chrMito 25318-25342      + |       Q0045
    ##   [8]  chrMito 26229-26701      + |       Q0045
    ## 
    ## $Q0050 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames      ranges strand |  Name
    ##   [1]  chrMito 13818-16322      + | Q0050
    ## 
    ## $Q0055 
    ## GRanges object with 2 ranges and 1 metadata column:
    ##       seqnames      ranges strand |  Name
    ##   [1]  chrMito 13818-13986      + | Q0055
    ##   [2]  chrMito 16435-18830      + | Q0055
    ## 
    ## ...
    ## <4893 more elements>
    ## -------
    ## seqinfo: 18 sequences from an unspecified genome; no seqlengths

Get TE from the the paper
=========================

``` r
te_paper <- read_tsv("../annotations/GSE53313_Cerevisiae_RNA_RPF.txt", skip = 2) %>% 
  setNames(c("id", "total", "ribo", "te")) %>% 
  filter(!is.na(te) & ribo > 5 & total > 5) %>%
  print()
```

    ## # A tibble: 4,442 x 4
    ##    id       total    ribo      te
    ##    <chr>    <dbl>   <dbl>   <dbl>
    ##  1 YAL001C   34.5   18.0  -0.942 
    ##  2 YAL002W   19.8    5.37 -1.88  
    ##  3 YAL003W 4168.  6818.    0.710 
    ##  4 YAL007C  193.   195.    0.0135
    ##  5 YAL008W   44.7   85.7   0.941 
    ##  6 YAL009W   31.0   17.5  -0.829 
    ##  7 YAL010C   14.9   10.3  -0.536 
    ##  8 YAL011W   16.3   10.4  -0.652 
    ##  9 YAL012W  518.  1103.    1.09  
    ## 10 YAL013W   23.6   40.2   0.768 
    ## # ... with 4,432 more rows

Load RQC stalls for joining with high TE genes
==============================================

``` r
rqc_stalls <- read_tsv("../../rqc_stalls_in_yeast_orfs/tables/ngrams_annotated.tsv") %>% 
  mutate(stall = if_else(stall %in% c("KR", "P"), "KPR", stall)) %>%
  print()
```

    ## # A tibble: 1,251 x 7
    ##    id     ngram   ngram_weight   pos stall gene  Note                     
    ##    <chr>  <chr>          <int> <int> <chr> <chr> <chr>                    
    ##  1 YHR13… RRRRRR…           10   321 KPR   <NA>  Putative protein of unkn…
    ##  2 YIL15… PPPPPP…           10   768 KPR   BNR1  Formin, nucleates the fo…
    ##  3 YNL27… PPPPPP…           10  1238 KPR   BNI1  Formin, nucleates the fo…
    ##  4 YOR01… KKKKKK…           10   712 KPR   <NA>  Protein of unknown funct…
    ##  5 YBL09… KKKKNK…            9    40 KPR   MAP2  Methionine aminopeptidas…
    ##  6 YDL14… APPPPP…            9   473 KPR   LDB17 Protein involved in the …
    ##  7 YDL17… KEKKKK…            9   273 KPR   PAR32 Putative protein of unkn…
    ##  8 YDR17… KKKKKK…            9   227 KPR   HMO1  Chromatin associated hig…
    ##  9 YDR33… KKKRGR…            9   231 KPR   SWR1  Swi2/Snf2-related ATPase…
    ## 10 YHL01… KRKDKK…            9   174 KPR   APM2  Protein of unknown funct…
    ## # ... with 1,241 more rows

Find typical (median) location of stall for use in control
==========================================================

``` r
median_stall_loc <- median(rqc_stalls[["pos"]]) %>% 
  print()
```

    ## [1] 215

Convert RQC stalls to genomic coordinates
=========================================

``` r
rqc_stalls_coords <- rqc_stalls %>% 
  mutate(seqname = id, start = pos*3 + 1) %>%
  mutate(end = start) %>%
  select(seqname, start, end, stall, id) %>%
  GRanges() %>%
  mapFromTranscripts(tx) %>% 
  # get rid of mitochondrial sequence
  filter(seqnames != "chrMito") %>% 
  mutate(id = rqc_stalls[xHits, "id"], stall = rqc_stalls[xHits, "stall"]) %>%
  select(-xHits, -transcriptsHits)

# check that the mapping was done correctly
rqc_stalls_coords %>% 
  anchor_5p() %>% 
  stretch(29) %>% 
  getSeq(genome, .) %>% 
  translate()
```

    ##   A AAStringSet instance of length 1102
    ##        width seq
    ##    [1]    10 PPPPPPPPPP
    ##    [2]    10 PPPPPPPPPP
    ##    [3]    10 KKKKKKKKKK
    ##    [4]    10 KKKKNKKKKK
    ##    [5]    10 APPPPPPPPP
    ##    ...   ... ...
    ## [1098]    10 KSVKKRKIMK
    ## [1099]    10 KPKPTPPSPP
    ## [1100]    10 NTKKKSRAKK
    ## [1101]    10 KISNRLRKRR
    ## [1102]    10 KPRGRKGGRK

Load the Ribo-seq alignments
============================

We do not run the codecell below after the first time to save time.

``` r
aln <- readGAlignments("../processeddata/mono/accepted_hits.bam") %>% 
  print()
```

Trim the Ribo-seq alignments to the P-site and calculate coverage separately for + and - strands
================================================================================================

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

Load the RNA-seq alignments
===========================

We do not run the codecell below after the first time to save time.

``` r
aln <- readGAlignments("../processeddata/totalrz/accepted_hits.bam") %>% 
  print()
```

Trim the RNA-seq alignments to the P-site and calculate coverage separately for + and - strands
===============================================================================================

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

rtracklayer::export.bw(cvg_plus, "../processeddata/totalrz/cvg_plus.bw")
rtracklayer::export.bw(cvg_minus, "../processeddata/totalrz/cvg_minus.bw")
```

Load pre-computed coverage
==========================

``` r
ribo_cvg_plus <- rtracklayer::import.bw("../processeddata/mono/cvg_plus.bw") %>% 
  coverage(weight = "score")
ribo_cvg_minus <- rtracklayer::import.bw("../processeddata/mono/cvg_minus.bw") %>% 
  coverage(weight = "score")
total_cvg_plus <- rtracklayer::import.bw("../processeddata/totalrz/cvg_plus.bw") %>% 
  coverage(weight = "score")
total_cvg_minus <- rtracklayer::import.bw("../processeddata/totalrz/cvg_minus.bw") %>% 
  coverage(weight = "score")
```

``` r
cvg <- c('ribo' = c(GRanges(ribo_cvg_plus, strand = "+"), GRanges(ribo_cvg_minus, strand = "-")),
  'total' = c(GRanges(total_cvg_plus, strand = "+"), GRanges(total_cvg_minus, strand = "-"))
  ) %>% 
  GRangesList() %>% 
  as_tibble() %>% 
  select(-group) %>% 
  rename(sample = group_name) %>% 
  filter(score > 0) %>% 
  GRanges()
```

``` r
tx_counts <- mapToTranscripts(cvg, tx) %>% 
  GRanges() %>% 
  mutate(sample = cvg$sample[xHits], score = cvg$score[xHits]) %>% 
  as_tibble()
```

``` r
te_calculated <- tx_counts %>% 
  group_by(seqnames, sample) %>% 
  summarize(counts = sum(width * score)) %>% 
  ungroup() %>% 
  spread(sample, counts) %>% 
  filter(ribo > 100 & total > 100) %>%
  mutate(te = log2(ribo / total / sum(ribo) * sum(total))) %>%
  dplyr::rename(id = seqnames) %>% 
  print()
```

    ## # A tibble: 4,284 x 4
    ##    id        ribo total     te
    ##    <fct>    <dbl> <dbl>  <dbl>
    ##  1 YAL001C   3232   867 -1.18 
    ##  2 YAL002W   1122   559 -2.07 
    ##  3 YAL003W 223552 19006  0.478
    ##  4 YAL005C 701880 74785  0.152
    ##  5 YAL007C   6906   892 -0.126
    ##  6 YAL008W   2897   201  0.771
    ##  7 YAL009W    753   181 -1.02 
    ##  8 YAL010C    838   166 -0.743
    ##  9 YAL011W   1015   218 -0.859
    ## 10 YAL012W  66467  4470  0.816
    ## # ... with 4,274 more rows

Plot calculated TE vs. TE pre-calculated by Weinberg et al.
===========================================================

``` r
te_paper %>% 
  select(id, te) %>% 
  inner_join(te_calculated, by = "id") %>% 
  ggplot(aes(x = te.x, y = te.y)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(limits = c(-6.5, 2.5)) +
  scale_y_continuous(limits = c(-6.5, 2.5)) +
  geom_abline(intercept = 0, slope = 1)
```

![](plot_te_for_only_preceding_stall_region_files/figure-markdown_github/unnamed-chunk-10-1.png)

Look at high TE (paper) genes with potential RQC stalls
=======================================================

``` r
te_paper_stall_data <- te_paper %>% 
  left_join(rqc_stalls, by = "id") %>%
  arrange(desc(te)) %>%
  select(te, gene, pos, ngram, ngram_weight, everything()) %>%
  print()
```

    ## # A tibble: 4,442 x 10
    ##       te gene    pos ngram  ngram_weight id    total   ribo stall Note    
    ##    <dbl> <chr> <int> <chr>         <int> <chr> <dbl>  <dbl> <chr> <chr>   
    ##  1  1.98 <NA>     NA <NA>             NA YLR4…  80.8  318.  <NA>  <NA>    
    ##  2  1.84 <NA>     NA <NA>             NA YPR0… 363.  1296.  <NA>  <NA>    
    ##  3  1.76 <NA>     NA <NA>             NA YMR1… 913.  3084.  <NA>  <NA>    
    ##  4  1.71 <NA>     NA <NA>             NA YJL1…  62.0  203.  <NA>  <NA>    
    ##  5  1.70 <NA>     NA <NA>             NA YDL0…  68.8  224.  <NA>  <NA>    
    ##  6  1.69 NOP16   215 KRRLL…            7 YER0…  59.9  193.  KPR   Constit…
    ##  7  1.65 <NA>     NA <NA>             NA YKL0…  19.1   60.0 <NA>  <NA>    
    ##  8  1.63 <NA>     NA <NA>             NA YBR0… 123.   382.  <NA>  <NA>    
    ##  9  1.60 <NA>     NA <NA>             NA YKL1… 335.  1016.  <NA>  <NA>    
    ## 10  1.58 <NA>     NA <NA>             NA YOR2…  30.2   90.5 <NA>  <NA>    
    ## # ... with 4,432 more rows

Plot TE from paper as a function of stall strength
==================================================

``` r
plot_data <- te_paper_stall_data %>% 
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

![](plot_te_for_only_preceding_stall_region_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
ggsave("../figures/distribution_of_te_paper_for_rqc_stall_containing_saccer_genes.pdf")
```

Test if stall-containing genes have lower or higher TE from paper
=================================================================

``` r
wilcox.test(te ~ ngram_weight, data = plot_data, alternative = "two.sided") %>% 
  tidy() %>% 
  gather() %>% 
  knitr::kable()
```

| key         | value                                             |
|:------------|:--------------------------------------------------|
| statistic   | 2031666.5                                         |
| p.value     | 2.2570475508283e-17                               |
| method      | Wilcoxon rank sum test with continuity correction |
| alternative | two.sided                                         |

Plot difference in TE between stall-containing and remaining genes
==================================================================

``` r
plot_data %>% 
  group_by(ngram_weight) %>% 
  summarize(median(te)) %>% 
  knitr::kable()
```

| ngram\_weight |  median(te)|
|:--------------|-----------:|
| No stall      |  -0.1943524|
| Stall         |  -0.4444032|

Look at high TE (calculated) genes with potential RQC stalls
============================================================

``` r
te_calc_stall_data <- te_calculated %>% 
  left_join(rqc_stalls, by = "id") %>%
  arrange(desc(te)) %>%
  select(te, gene, pos, ngram, ngram_weight, everything()) %>%
  print()
```

    ## # A tibble: 4,284 x 10
    ##       te gene    pos ngram  ngram_weight id     ribo total stall Note     
    ##    <dbl> <chr> <int> <chr>         <int> <chr> <dbl> <dbl> <chr> <chr>    
    ##  1  2.02 <NA>     NA <NA>             NA YDL0… 66683  1944 <NA>  <NA>     
    ##  2  1.90 <NA>     NA <NA>             NA YLR4…  4606   146 <NA>  <NA>     
    ##  3  1.65 <NA>     NA <NA>             NA YPR0… 12927   489 <NA>  <NA>     
    ##  4  1.58 <NA>     NA <NA>             NA YDR4…  8761   348 <NA>  <NA>     
    ##  5  1.53 NOP16   215 KRRLL…            7 YER0…  7432   305 KPR   Constitu…
    ##  6  1.51 <NA>     NA <NA>             NA YDR0… 89386  3719 <NA>  <NA>     
    ##  7  1.50 <NA>     NA <NA>             NA YJL1…  6912   289 <NA>  <NA>     
    ##  8  1.45 <NA>     NA <NA>             NA YDL0…  3557   154 <NA>  <NA>     
    ##  9  1.38 <NA>     NA <NA>             NA YLR3… 71326  3240 <NA>  <NA>     
    ## 10  1.38 <NA>     NA <NA>             NA YDR3…  3033   138 <NA>  <NA>     
    ## # ... with 4,274 more rows

Plot calculated TE as a function of stall strength
==================================================

``` r
plot_data <- te_calc_stall_data %>% 
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

![](plot_te_for_only_preceding_stall_region_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
ggsave("../figures/distribution_of_te_calculated_for_rqc_stall_containing_saccer_genes.pdf")
```

Test if stall-containing genes have lower or higher calculated TE
=================================================================

``` r
wilcox.test(te ~ ngram_weight, data = plot_data, alternative = "two.sided") %>% 
  tidy() %>% 
  gather() %>% 
  knitr::kable()
```

| key         | value                                             |
|:------------|:--------------------------------------------------|
| statistic   | 1956806                                           |
| p.value     | 4.3544146969543e-17                               |
| method      | Wilcoxon rank sum test with continuity correction |
| alternative | two.sided                                         |

Plot difference in TE between stall-containing and remaining genes
==================================================================

``` r
plot_data %>% 
  group_by(ngram_weight) %>% 
  summarize(median(te)) %>% 
  knitr::kable()
```

| ngram\_weight |  median(te)|
|:--------------|-----------:|
| No stall      |  -0.3820747|
| Stall         |  -0.6813425|

Calculate TE for only the region preceding each stall
=====================================================

``` r
te_preceding_stalls <- tx_counts %>% 
  left_join(rqc_stalls %>% select(id, pos), by = c("seqnames" = "id")) %>% 
  mutate(pos = if_else(is.na(pos), median_stall_loc, pos)) %>% 
  filter(start < 3*pos + 1) %>% 
  group_by(seqnames, sample) %>% 
  summarize(counts = sum(width * score)) %>% 
  ungroup() %>% 
  spread(sample, counts) %>% 
  filter(ribo > 100 & total > 100) %>%
  mutate(te = log2(ribo / total / sum(ribo) * sum(total))) %>%
  dplyr::rename(id = seqnames) %>% 
  print()
```

    ## # A tibble: 3,622 x 4
    ##    id        ribo total     te
    ##    <chr>    <dbl> <dbl>  <dbl>
    ##  1 YAL001C    869   203 -1.19 
    ##  2 YAL003W 223552 19006  0.266
    ##  3 YAL005C 208604 19538  0.126
    ##  4 YAL007C   6906   891 -0.336
    ##  5 YAL008W   2897   201  0.559
    ##  6 YAL009W    595   171 -1.49 
    ##  7 YAL012W  32012  2135  0.616
    ##  8 YAL014C    986   270 -1.42 
    ##  9 YAL016W   4722   547 -0.180
    ## 10 YAL017W   1304   266 -0.996
    ## # ... with 3,612 more rows

Look at high TE (calculated only preceding stall) genes with potential RQC stalls
=================================================================================

``` r
te_preceding_stall_data <- te_preceding_stalls %>% 
  left_join(rqc_stalls, by = "id") %>%
  arrange(desc(te)) %>%
  select(te, gene, pos, ngram, ngram_weight, everything()) %>%
  print()
```

    ## # A tibble: 3,622 x 10
    ##       te gene    pos ngram  ngram_weight id     ribo total stall Note     
    ##    <dbl> <chr> <int> <chr>         <int> <chr> <dbl> <dbl> <chr> <chr>    
    ##  1  1.81 <NA>     NA <NA>             NA YDL0… 66683  1944 <NA>  <NA>     
    ##  2  1.70 HTB1     28 DGKKR…            6 YDR2… 10546   331 KPR   Histone …
    ##  3  1.69 <NA>     NA <NA>             NA YLR4…  4606   146 <NA>  <NA>     
    ##  4  1.49 <NA>     NA <NA>             NA YDR0… 60140  2195 <NA>  <NA>     
    ##  5  1.44 NOP16   215 KRRLL…            7 YER0…  6866   259 KPR   Constitu…
    ##  6  1.43 <NA>     NA <NA>             NA YPR0… 12927   489 <NA>  <NA>     
    ##  7  1.36 <NA>     NA <NA>             NA YDR4…  8761   348 <NA>  <NA>     
    ##  8  1.29 <NA>     NA <NA>             NA YJL1…  6912   289 <NA>  <NA>     
    ##  9  1.24 <NA>     NA <NA>             NA YDL0…  3557   154 <NA>  <NA>     
    ## 10  1.21 NPL3     94 PQPYY…            7 YDR4…  8231   363 KPR   RNA-bind…
    ## # ... with 3,612 more rows

Plot calculated TE for region preceding stalls as a function of stall strength
==============================================================================

``` r
plot_data <- te_preceding_stall_data %>% 
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

![](plot_te_for_only_preceding_stall_region_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
ggsave("../figures/distribution_of_te_preceding_stalls_for_rqc_stall_containing_saccer_genes.pdf")
```

Test if stall-containing genes have lower or higher calculated TE preceding stall
=================================================================================

``` r
wilcox.test(te ~ ngram_weight, data = plot_data, alternative = "two.sided") %>% 
  tidy() %>% 
  gather() %>% 
  knitr::kable()
```

| key         | value                                             |
|:------------|:--------------------------------------------------|
| statistic   | 1362395.5                                         |
| p.value     | 2.10316183148019e-26                              |
| method      | Wilcoxon rank sum test with continuity correction |
| alternative | two.sided                                         |

Plot difference in TE between stall-containing and remaining genes
==================================================================

``` r
plot_data %>% 
  group_by(ngram_weight) %>% 
  summarize(median(te)) %>% 
  knitr::kable()
```

| ngram\_weight |  median(te)|
|:--------------|-----------:|
| No stall      |  -0.3310395|
| Stall         |  -0.6506406|

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
    ##  [3] biobroom_1.12.1             broom_0.5.0                
    ##  [5] plyranges_1.0.3             forcats_0.3.0              
    ##  [7] stringr_1.3.1               dplyr_0.7.6                
    ##  [9] purrr_0.2.5                 readr_1.1.1                
    ## [11] tidyr_0.8.1                 tibble_1.4.2               
    ## [13] ggplot2_3.0.0               tidyverse_1.2.1            
    ## [15] GenomicFeatures_1.32.0      AnnotationDbi_1.42.1       
    ## [17] GenomicAlignments_1.16.0    Rsamtools_1.32.2           
    ## [19] Biostrings_2.48.0           XVector_0.20.0             
    ## [21] SummarizedExperiment_1.10.1 DelayedArray_0.6.2         
    ## [23] BiocParallel_1.14.2         matrixStats_0.54.0         
    ## [25] Biobase_2.40.0              GenomicRanges_1.32.6       
    ## [27] GenomeInfoDb_1.16.0         IRanges_2.14.10            
    ## [29] S4Vectors_0.18.3            BiocGenerics_0.26.0        
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
    ## [31] highr_0.7                              
    ## [32] BSgenome_1.48.0                        
    ## [33] rlang_0.2.1                            
    ## [34] readxl_1.1.0                           
    ## [35] rstudioapi_0.7                         
    ## [36] RSQLite_2.1.1                          
    ## [37] bindr_0.1.1                            
    ## [38] jsonlite_1.5                           
    ## [39] RCurl_1.95-4.11                        
    ## [40] magrittr_1.5                           
    ## [41] GenomeInfoDbData_1.1.0                 
    ## [42] Matrix_1.2-14                          
    ## [43] Rcpp_0.12.18                           
    ## [44] munsell_0.5.0                          
    ## [45] fansi_0.2.3                            
    ## [46] stringi_1.2.4                          
    ## [47] yaml_2.2.0                             
    ## [48] zlibbioc_1.26.0                        
    ## [49] plyr_1.8.4                             
    ## [50] grid_3.5.1                             
    ## [51] blob_1.1.1                             
    ## [52] crayon_1.3.4                           
    ## [53] lattice_0.20-35                        
    ## [54] haven_1.1.2                            
    ## [55] hms_0.4.2                              
    ## [56] knitr_1.20                             
    ## [57] pillar_1.3.0                           
    ## [58] biomaRt_2.36.1                         
    ## [59] XML_3.98-1.12                          
    ## [60] glue_1.3.0                             
    ## [61] evaluate_0.11                          
    ## [62] modelr_0.1.2                           
    ## [63] cellranger_1.1.0                       
    ## [64] gtable_0.2.0                           
    ## [65] assertthat_0.2.0                       
    ## [66] memoise_1.1.0

Source data for S5 Fig panel B
==============================

``` r
plot_data %>% 
  rename(x = ngram_weight, y = te) %>% 
  select(id, x, y) %>% 
  mutate_if(is.numeric, funs(round(., 3))) %>% 
  knitr::kable()
```

| id        | x        |       y|
|:----------|:---------|-------:|
| YDL061C   | No stall |   1.810|
| YDR224C   | Stall    |   1.704|
| YLR438C-A | No stall |   1.690|
| YDR077W   | No stall |   1.486|
| YER002W   | Stall    |   1.439|
| YPR036W-A | No stall |   1.434|
| YDR461W   | No stall |   1.364|
| YJL143W   | No stall |   1.290|
| YDL012C   | No stall |   1.240|
| YDR432W   | Stall    |   1.213|
| YKL164C   | No stall |   1.198|
| YLR388W   | No stall |   1.170|
| YDR322C-A | No stall |   1.168|
| YNL215W   | Stall    |   1.147|
| YBR016W   | No stall |   1.118|
| YJL159W   | No stall |   1.099|
| YKL117W   | No stall |   1.098|
| YNL145W   | No stall |   1.076|
| YHL034C   | No stall |   1.074|
| YHR001W-A | No stall |   1.070|
| YOR257W   | No stall |   1.048|
| YJL158C   | No stall |   1.043|
| YEL017C-A | No stall |   1.015|
| YMR256C   | No stall |   0.998|
| YCR024C-A | No stall |   0.992|
| YLR390W-A | No stall |   0.971|
| YGR008C   | No stall |   0.957|
| YLR025W   | No stall |   0.954|
| YDL014W   | No stall |   0.949|
| YDL067C   | No stall |   0.941|
| YER006W   | Stall    |   0.941|
| YIL008W   | No stall |   0.929|
| YBR082C   | No stall |   0.920|
| YPR062W   | No stall |   0.909|
| YLR038C   | No stall |   0.901|
| YBL092W   | Stall    |   0.896|
| YJL166W   | No stall |   0.892|
| YGR280C   | Stall    |   0.887|
| YFL010C   | No stall |   0.886|
| YJL115W   | No stall |   0.885|
| YJR048W   | No stall |   0.875|
| YDL130W-A | No stall |   0.868|
| YCL035C   | No stall |   0.867|
| YDL184C   | Stall    |   0.866|
| YNL070W   | No stall |   0.862|
| YHR053C   | No stall |   0.861|
| YHR055C   | No stall |   0.861|
| YLR435W   | Stall    |   0.857|
| YML058W   | No stall |   0.849|
| YBR173C   | No stall |   0.846|
| YBR009C   | No stall |   0.836|
| YPR020W   | No stall |   0.826|
| YNL007C   | Stall    |   0.823|
| YJL079C   | No stall |   0.817|
| YCR028C-A | No stall |   0.810|
| YCL050C   | No stall |   0.800|
| YML129C   | No stall |   0.798|
| YAL025C   | Stall    |   0.797|
| YHR089C   | Stall    |   0.794|
| YOR089C   | No stall |   0.789|
| YLR147C   | No stall |   0.787|
| YLL009C   | No stall |   0.777|
| YJR014W   | Stall    |   0.776|
| YIR037W   | No stall |   0.776|
| YLR395C   | No stall |   0.772|
| YHR143W-A | No stall |   0.771|
| YKL096W-A | No stall |   0.771|
| YGR037C   | No stall |   0.767|
| YBR067C   | No stall |   0.767|
| YJR047C   | No stall |   0.763|
| YKL192C   | No stall |   0.763|
| YOR252W   | Stall    |   0.752|
| YNL024C-A | No stall |   0.751|
| YGL029W   | Stall    |   0.751|
| YNL079C   | No stall |   0.750|
| YER100W   | No stall |   0.749|
| YDR162C   | No stall |   0.734|
| YOL109W   | No stall |   0.730|
| YJR125C   | Stall    |   0.729|
| YPL249C-A | No stall |   0.728|
| YOR265W   | No stall |   0.725|
| YPL170W   | No stall |   0.722|
| YDR412W   | Stall    |   0.716|
| YMR039C   | Stall    |   0.714|
| YFR049W   | No stall |   0.711|
| YOL077W-A | No stall |   0.707|
| YCL037C   | No stall |   0.706|
| YER011W   | No stall |   0.705|
| YMR173W   | No stall |   0.704|
| YHR193C   | No stall |   0.700|
| YNL030W   | No stall |   0.689|
| YJR069C   | No stall |   0.688|
| YHR072W-A | No stall |   0.687|
| YOR159C   | No stall |   0.679|
| YOL155C   | No stall |   0.677|
| YDR299W   | No stall |   0.671|
| YEL034W   | No stall |   0.669|
| YFR032C-A | No stall |   0.668|
| YNL111C   | No stall |   0.667|
| YLR009W   | Stall    |   0.667|
| YIR034C   | No stall |   0.661|
| YKR043C   | No stall |   0.659|
| YPL163C   | No stall |   0.656|
| YOR332W   | No stall |   0.654|
| YJR010C-A | No stall |   0.650|
| YNL055C   | No stall |   0.648|
| YOR197W   | Stall    |   0.647|
| YKL156W   | No stall |   0.645|
| YPL206C   | No stall |   0.644|
| YPR035W   | Stall    |   0.642|
| YOL005C   | No stall |   0.641|
| YLR109W   | No stall |   0.640|
| YNL061W   | Stall    |   0.636|
| YMR251W-A | No stall |   0.633|
| YGR254W   | No stall |   0.628|
| YJL123C   | No stall |   0.628|
| YPR183W   | No stall |   0.619|
| YMR002W   | No stall |   0.618|
| YPL061W   | No stall |   0.617|
| YAL012W   | No stall |   0.616|
| YLR264W   | No stall |   0.614|
| YDR032C   | No stall |   0.612|
| YLR200W   | No stall |   0.609|
| YGR191W   | No stall |   0.604|
| YDR468C   | No stall |   0.604|
| YNL015W   | No stall |   0.604|
| YOR327C   | No stall |   0.598|
| YMR203W   | No stall |   0.596|
| YKL138C-A | No stall |   0.594|
| YGL044C   | No stall |   0.590|
| YBL030C   | No stall |   0.586|
| YGR279C   | No stall |   0.581|
| YOR382W   | No stall |   0.579|
| YJR060W   | No stall |   0.575|
| YKL054C   | No stall |   0.574|
| YMR295C   | Stall    |   0.573|
| YML081C-A | No stall |   0.571|
| YKR014C   | No stall |   0.570|
| YEL020W-A | No stall |   0.569|
| YFR001W   | Stall    |   0.568|
| YJL189W   | Stall    |   0.567|
| YHR021C   | No stall |   0.567|
| YHR174W   | No stall |   0.564|
| YNL096C   | No stall |   0.563|
| YAL008W   | No stall |   0.559|
| YLR058C   | No stall |   0.558|
| YDR171W   | No stall |   0.558|
| YMR091C   | Stall    |   0.558|
| YDR510W   | No stall |   0.556|
| YNL052W   | No stall |   0.554|
| YDR493W   | No stall |   0.553|
| YMR072W   | No stall |   0.552|
| YDR384C   | No stall |   0.551|
| YOR052C   | Stall    |   0.548|
| YOR045W   | No stall |   0.545|
| YOR383C   | No stall |   0.540|
| YPL183W-A | No stall |   0.539|
| YPR125W   | No stall |   0.538|
| YBR109C   | No stall |   0.537|
| YDR353W   | No stall |   0.536|
| YJL179W   | No stall |   0.535|
| YPL111W   | No stall |   0.535|
| YHR152W   | Stall    |   0.533|
| YBR058C-A | No stall |   0.531|
| YPL037C   | No stall |   0.529|
| YER056C-A | Stall    |   0.527|
| YGL106W   | No stall |   0.526|
| YDL075W   | Stall    |   0.525|
| YNL147W   | No stall |   0.524|
| YIL051C   | No stall |   0.523|
| YBR046C   | No stall |   0.523|
| YNL166C   | No stall |   0.523|
| YMR311C   | No stall |   0.521|
| YLR170C   | No stall |   0.519|
| YKL096W   | No stall |   0.517|
| YDR276C   | No stall |   0.516|
| YKL163W   | No stall |   0.514|
| YJR104C   | No stall |   0.509|
| YDR158W   | No stall |   0.506|
| YDR363W-A | No stall |   0.503|
| YDR092W   | No stall |   0.503|
| YBL028C   | No stall |   0.500|
| YDL166C   | No stall |   0.499|
| YLR355C   | No stall |   0.497|
| YOR027W   | No stall |   0.496|
| YLR061W   | No stall |   0.495|
| YBL091C   | Stall    |   0.495|
| YNR046W   | No stall |   0.494|
| YJR045C   | No stall |   0.491|
| YCL028W   | No stall |   0.490|
| YNL255C   | No stall |   0.490|
| YLR325C   | No stall |   0.490|
| YLR350W   | No stall |   0.490|
| YKL067W   | No stall |   0.489|
| YLR185W   | No stall |   0.488|
| YDR513W   | No stall |   0.487|
| YBR010W   | No stall |   0.485|
| YER048C   | No stall |   0.484|
| YDL022W   | No stall |   0.484|
| YDR178W   | No stall |   0.482|
| YLR110C   | No stall |   0.481|
| YHL021C   | No stall |   0.481|
| YDL153C   | Stall    |   0.480|
| YNL160W   | No stall |   0.478|
| YFR011C   | No stall |   0.475|
| YMR194W   | Stall    |   0.474|
| YGL070C   | No stall |   0.473|
| YHR059W   | No stall |   0.472|
| YML028W   | No stall |   0.472|
| YHR039C-A | No stall |   0.472|
| YKR025W   | No stall |   0.471|
| YHR081W   | Stall    |   0.471|
| YDR177W   | No stall |   0.468|
| YKR013W   | No stall |   0.467|
| YGR019W   | No stall |   0.467|
| YDR113C   | No stall |   0.458|
| YPL098C   | No stall |   0.458|
| YIL033C   | No stall |   0.455|
| YAL035W   | Stall    |   0.455|
| YJL140W   | Stall    |   0.454|
| YER120W   | No stall |   0.454|
| YER143W   | No stall |   0.453|
| YDR272W   | No stall |   0.453|
| YHL015W   | No stall |   0.451|
| YNL036W   | No stall |   0.448|
| YNL031C   | No stall |   0.445|
| YHR087W   | No stall |   0.445|
| YGR105W   | No stall |   0.445|
| YBR126C   | No stall |   0.444|
| YPL127C   | No stall |   0.444|
| YLR406C   | Stall    |   0.443|
| YBR072W   | No stall |   0.443|
| YKR001C   | No stall |   0.442|
| YNL064C   | No stall |   0.439|
| YLR055C   | No stall |   0.438|
| YER025W   | No stall |   0.438|
| YDR462W   | No stall |   0.438|
| YHR051W   | No stall |   0.433|
| YFR053C   | No stall |   0.432|
| YFL017W-A | No stall |   0.432|
| YDL065C   | No stall |   0.431|
| YGL187C   | No stall |   0.431|
| YDR167W   | No stall |   0.430|
| YDR155C   | No stall |   0.430|
| YOR298C-A | No stall |   0.429|
| YER003C   | No stall |   0.427|
| YLL039C   | No stall |   0.426|
| YOL111C   | No stall |   0.425|
| YHR018C   | No stall |   0.424|
| YLR295C   | No stall |   0.423|
| YNR054C   | No stall |   0.422|
| YDL110C   | No stall |   0.422|
| YNL208W   | No stall |   0.421|
| YLR438W   | No stall |   0.421|
| YGR209C   | No stall |   0.420|
| YGL172W   | No stall |   0.420|
| YCL064C   | No stall |   0.419|
| YMR193W   | Stall    |   0.417|
| YNL274C   | No stall |   0.417|
| YOL032W   | No stall |   0.417|
| YNL175C   | Stall    |   0.416|
| YPL117C   | No stall |   0.416|
| YGR167W   | No stall |   0.415|
| YML074C   | Stall    |   0.415|
| YPL271W   | No stall |   0.415|
| YIL052C   | No stall |   0.414|
| YIL123W   | No stall |   0.412|
| YNL110C   | No stall |   0.412|
| YDR172W   | No stall |   0.410|
| YHR179W   | No stall |   0.410|
| YGL191W   | No stall |   0.409|
| YDL135C   | No stall |   0.408|
| YOR094W   | No stall |   0.406|
| YGR063C   | No stall |   0.406|
| YCR047C   | Stall    |   0.403|
| YGR183C   | No stall |   0.403|
| YER074W   | Stall    |   0.403|
| YKL035W   | No stall |   0.398|
| YGL030W   | No stall |   0.398|
| YLR150W   | No stall |   0.397|
| YHR092C   | No stall |   0.397|
| YGL244W   | Stall    |   0.397|
| YDR002W   | No stall |   0.393|
| YPL267W   | No stall |   0.392|
| YJR034W   | No stall |   0.392|
| YBR249C   | No stall |   0.392|
| YCR088W   | Stall    |   0.391|
| YKL009W   | No stall |   0.389|
| YKR081C   | No stall |   0.389|
| YLR208W   | No stall |   0.389|
| YNL135C   | No stall |   0.387|
| YGR152C   | No stall |   0.387|
| YKR035W-A | No stall |   0.386|
| YLR029C   | Stall    |   0.384|
| YDL161W   | Stall    |   0.382|
| YDR225W   | No stall |   0.381|
| YOR317W   | No stall |   0.379|
| YGR215W   | Stall    |   0.378|
| YOL086C   | No stall |   0.378|
| YIL062C   | No stall |   0.376|
| YPL225W   | No stall |   0.376|
| YPL091W   | No stall |   0.375|
| YDR368W   | No stall |   0.370|
| YDR071C   | No stall |   0.369|
| YDL043C   | No stall |   0.368|
| YOL143C   | No stall |   0.367|
| YGL226C-A | No stall |   0.367|
| YDR129C   | No stall |   0.367|
| YDR267C   | No stall |   0.367|
| YJL190C   | No stall |   0.366|
| YBR290W   | No stall |   0.366|
| YDL168W   | No stall |   0.366|
| YDR512C   | No stall |   0.365|
| YGR275W   | No stall |   0.364|
| YKR042W   | No stall |   0.363|
| YML078W   | No stall |   0.360|
| YML063W   | Stall    |   0.359|
| YER112W   | No stall |   0.358|
| YOR293W   | No stall |   0.355|
| YJL136C   | No stall |   0.355|
| YDR055W   | No stall |   0.354|
| YLR262C   | No stall |   0.354|
| YER031C   | No stall |   0.353|
| YGL026C   | No stall |   0.351|
| YHL031C   | No stall |   0.350|
| YDR500C   | No stall |   0.349|
| YBR154C   | No stall |   0.348|
| YLR367W   | No stall |   0.348|
| YBR078W   | No stall |   0.348|
| YEL037C   | No stall |   0.347|
| YDR156W   | No stall |   0.346|
| YKL087C   | No stall |   0.346|
| YER057C   | No stall |   0.345|
| YOR198C   | Stall    |   0.344|
| YGR076C   | No stall |   0.344|
| YKL172W   | Stall    |   0.344|
| YNL004W   | No stall |   0.343|
| YDR381W   | Stall    |   0.343|
| YGL037C   | No stall |   0.338|
| YKL216W   | No stall |   0.337|
| YBR011C   | No stall |   0.335|
| YDR424C   | No stall |   0.335|
| YLR448W   | Stall    |   0.333|
| YML094W   | No stall |   0.330|
| YIL069C   | Stall    |   0.328|
| YDR190C   | No stall |   0.328|
| YMR276W   | No stall |   0.328|
| YDR099W   | No stall |   0.327|
| YKL032C   | Stall    |   0.327|
| YMR121C   | Stall    |   0.326|
| YGR080W   | Stall    |   0.325|
| YHR146W   | Stall    |   0.324|
| YDR098C   | No stall |   0.323|
| YJL026W   | No stall |   0.323|
| YGR180C   | No stall |   0.323|
| YCR004C   | No stall |   0.323|
| YLR153C   | No stall |   0.321|
| YMR230W   | No stall |   0.321|
| YPL129W   | No stall |   0.320|
| YBL050W   | No stall |   0.320|
| YIL070C   | No stall |   0.319|
| YOR239W   | Stall    |   0.319|
| YLR287C-A | Stall    |   0.314|
| YLR270W   | No stall |   0.314|
| YGR251W   | Stall    |   0.313|
| YLL050C   | No stall |   0.313|
| YPL198W   | No stall |   0.310|
| YNL113W   | No stall |   0.310|
| YIL106W   | No stall |   0.310|
| YML073C   | Stall    |   0.310|
| YBL003C   | No stall |   0.309|
| YNL131W   | No stall |   0.309|
| YNR019W   | No stall |   0.309|
| YOL040C   | No stall |   0.305|
| YFL034C-A | No stall |   0.304|
| YGR118W   | No stall |   0.303|
| YPR043W   | No stall |   0.302|
| YDR064W   | No stall |   0.302|
| YAL044C   | No stall |   0.302|
| YLR293C   | No stall |   0.301|
| YBR035C   | No stall |   0.301|
| YFL005W   | No stall |   0.300|
| YPL239W   | No stall |   0.299|
| YDR328C   | No stall |   0.299|
| YLR051C   | Stall    |   0.298|
| YJR070C   | No stall |   0.297|
| YOR167C   | No stall |   0.296|
| YOR103C   | No stall |   0.296|
| YER087C-B | No stall |   0.295|
| YKL085W   | No stall |   0.295|
| YKL142W   | No stall |   0.295|
| YNL312W   | No stall |   0.293|
| YER090W   | No stall |   0.292|
| YBR137W   | No stall |   0.292|
| YPR041W   | No stall |   0.290|
| YLL028W   | No stall |   0.288|
| YEL046C   | No stall |   0.285|
| YDR289C   | No stall |   0.285|
| YMR263W   | Stall    |   0.284|
| YKR094C   | No stall |   0.282|
| YML001W   | No stall |   0.282|
| YBR164C   | No stall |   0.282|
| YMR071C   | No stall |   0.282|
| YDR296W   | Stall    |   0.282|
| YER072W   | No stall |   0.281|
| YJL062W-A | No stall |   0.280|
| YOL058W   | No stall |   0.280|
| YOR210W   | No stall |   0.277|
| YOR142W   | No stall |   0.277|
| YBL036C   | No stall |   0.276|
| YBR088C   | No stall |   0.276|
| YCL005W-A | No stall |   0.276|
| YBR247C   | No stall |   0.275|
| YOL121C   | No stall |   0.274|
| YJR105W   | No stall |   0.274|
| YLR133W   | No stall |   0.273|
| YOR294W   | Stall    |   0.270|
| YDR037W   | No stall |   0.268|
| YAL003W   | No stall |   0.266|
| YLR244C   | Stall    |   0.266|
| YLR180W   | No stall |   0.266|
| YLR321C   | Stall    |   0.265|
| YDR012W   | No stall |   0.264|
| YDR342C   | No stall |   0.264|
| YDR343C   | No stall |   0.264|
| YCR053W   | No stall |   0.263|
| YDR045C   | No stall |   0.262|
| YLR449W   | Stall    |   0.261|
| YMR286W   | No stall |   0.261|
| YPL013C   | No stall |   0.261|
| YMR318C   | No stall |   0.261|
| YBR031W   | No stall |   0.261|
| YPL004C   | No stall |   0.260|
| YIL148W   | No stall |   0.260|
| YBR135W   | No stall |   0.260|
| YNL037C   | No stall |   0.259|
| YLR370C   | No stall |   0.259|
| YOR006C   | No stall |   0.258|
| YHR132C   | No stall |   0.258|
| YKL042W   | Stall    |   0.258|
| YKL137W   | No stall |   0.257|
| YDR502C   | No stall |   0.257|
| YGR005C   | Stall    |   0.256|
| YFR004W   | No stall |   0.256|
| YML064C   | No stall |   0.254|
| YOR074C   | No stall |   0.254|
| YPR182W   | No stall |   0.253|
| YGR027C   | Stall    |   0.253|
| YGL076C   | No stall |   0.252|
| YML008C   | No stall |   0.251|
| YLR421C   | No stall |   0.250|
| YKR006C   | Stall    |   0.250|
| YDR050C   | No stall |   0.249|
| YJR025C   | No stall |   0.245|
| YLR418C   | No stall |   0.245|
| YBR251W   | No stall |   0.245|
| YMR315W   | No stall |   0.245|
| YGL222C   | No stall |   0.245|
| YPL235W   | No stall |   0.244|
| YER050C   | No stall |   0.243|
| YIL020C   | No stall |   0.242|
| YER063W   | No stall |   0.241|
| YDL160C   | No stall |   0.240|
| YGR181W   | No stall |   0.237|
| YLR437C   | No stall |   0.237|
| YBL069W   | No stall |   0.237|
| YER019C-A | No stall |   0.237|
| YNL301C   | No stall |   0.236|
| YER159C   | No stall |   0.236|
| YBL015W   | No stall |   0.235|
| YCR073W-A | Stall    |   0.235|
| YPR108W   | No stall |   0.234|
| YPL131W   | No stall |   0.234|
| YML126C   | No stall |   0.233|
| YDR148C   | No stall |   0.233|
| YGR141W   | No stall |   0.233|
| YPR166C   | No stall |   0.232|
| YIL138C   | No stall |   0.232|
| YBR248C   | No stall |   0.231|
| YOR020C   | No stall |   0.230|
| YJR135W-A | No stall |   0.230|
| YPL143W   | No stall |   0.230|
| YER117W   | Stall    |   0.230|
| YJL167W   | No stall |   0.229|
| YBR026C   | No stall |   0.228|
| YIL053W   | No stall |   0.226|
| YBR230C   | No stall |   0.226|
| YDL236W   | No stall |   0.225|
| YHR183W   | No stall |   0.223|
| YJL203W   | Stall    |   0.223|
| YOR259C   | No stall |   0.223|
| YOL038W   | No stall |   0.223|
| YDL125C   | No stall |   0.222|
| YDR086C   | No stall |   0.220|
| YBR244W   | No stall |   0.219|
| YBL071W-A | No stall |   0.219|
| YDL045W-A | No stall |   0.219|
| YDR481C   | Stall    |   0.218|
| YEL024W   | No stall |   0.218|
| YPR163C   | No stall |   0.217|
| YJL148W   | Stall    |   0.215|
| YDR044W   | No stall |   0.215|
| YKL150W   | No stall |   0.215|
| YBR014C   | No stall |   0.215|
| YAL038W   | No stall |   0.214|
| YLL026W   | No stall |   0.212|
| YDR023W   | No stall |   0.211|
| YNR044W   | No stall |   0.211|
| YPL059W   | No stall |   0.210|
| YER042W   | No stall |   0.210|
| YHR148W   | No stall |   0.210|
| YGR192C   | No stall |   0.209|
| YPL250C   | No stall |   0.207|
| YHR008C   | No stall |   0.206|
| YDL208W   | No stall |   0.205|
| YOR374W   | No stall |   0.205|
| YLL023C   | No stall |   0.204|
| YBR061C   | No stall |   0.203|
| YBL033C   | No stall |   0.203|
| YLR167W   | Stall    |   0.202|
| YDR226W   | No stall |   0.202|
| YMR143W   | No stall |   0.202|
| YPR132W   | No stall |   0.201|
| YML105C   | No stall |   0.201|
| YOR136W   | No stall |   0.201|
| YDR377W   | No stall |   0.200|
| YDR074W   | No stall |   0.200|
| YDR392W   | No stall |   0.200|
| YOR234C   | No stall |   0.199|
| YLL045C   | No stall |   0.199|
| YOR189W   | No stall |   0.199|
| YDR309C   | No stall |   0.197|
| YDL004W   | No stall |   0.197|
| YML069W   | No stall |   0.197|
| YPR065W   | Stall    |   0.196|
| YGR130C   | No stall |   0.194|
| YGL031C   | Stall    |   0.193|
| YGR285C   | Stall    |   0.193|
| YDL078C   | No stall |   0.193|
| YGR085C   | Stall    |   0.193|
| YER009W   | No stall |   0.192|
| YML070W   | No stall |   0.192|
| YNR043W   | No stall |   0.192|
| YNL149C   | Stall    |   0.192|
| YMR226C   | No stall |   0.191|
| YGR234W   | No stall |   0.190|
| YJL068C   | No stall |   0.190|
| YKL122C   | Stall    |   0.189|
| YIL118W   | No stall |   0.189|
| YNL151C   | No stall |   0.188|
| YIR001C   | No stall |   0.188|
| YER062C   | No stall |   0.187|
| YLR179C   | No stall |   0.187|
| YER074W-A | No stall |   0.187|
| YIL098C   | No stall |   0.186|
| YBL041W   | No stall |   0.186|
| YOR375C   | No stall |   0.186|
| YMR038C   | No stall |   0.185|
| YNL231C   | No stall |   0.183|
| YGR078C   | No stall |   0.183|
| YGR136W   | No stall |   0.181|
| YER052C   | No stall |   0.181|
| YFR050C   | No stall |   0.180|
| YBL045C   | No stall |   0.179|
| YLR043C   | No stall |   0.179|
| YOR287C   | Stall    |   0.179|
| YGR193C   | No stall |   0.177|
| YKL103C   | No stall |   0.177|
| YPL135W   | No stall |   0.176|
| YBL021C   | No stall |   0.175|
| YBR252W   | No stall |   0.175|
| YMR008C   | No stall |   0.175|
| YGL043W   | No stall |   0.174|
| YJR139C   | No stall |   0.174|
| YDR492W   | No stall |   0.173|
| YDL082W   | Stall    |   0.173|
| YFL039C   | No stall |   0.172|
| YIL019W   | Stall    |   0.172|
| YNL162W   | No stall |   0.171|
| YJR123W   | No stall |   0.171|
| YNL302C   | No stall |   0.170|
| YML092C   | No stall |   0.169|
| YLR196W   | Stall    |   0.169|
| YEL044W   | No stall |   0.169|
| YLR026C   | No stall |   0.169|
| YIL154C   | Stall    |   0.169|
| YHR190W   | No stall |   0.169|
| YGR187C   | No stall |   0.168|
| YFL022C   | No stall |   0.166|
| YML009C   | Stall    |   0.166|
| YGL202W   | No stall |   0.166|
| YNL067W   | No stall |   0.166|
| YKR057W   | No stall |   0.165|
| YHR038W   | No stall |   0.165|
| YMR112C   | No stall |   0.164|
| YHR020W   | No stall |   0.163|
| YOL039W   | No stall |   0.163|
| YMR298W   | No stall |   0.161|
| YDL083C   | No stall |   0.161|
| YOR247W   | No stall |   0.161|
| YOR063W   | No stall |   0.161|
| YIR012W   | No stall |   0.161|
| YER007C-A | No stall |   0.161|
| YOL144W   | Stall    |   0.160|
| YOL151W   | No stall |   0.160|
| YMR044W   | Stall    |   0.159|
| YHR094C   | No stall |   0.159|
| YLR439W   | No stall |   0.158|
| YAR015W   | No stall |   0.158|
| YJL001W   | No stall |   0.157|
| YER177W   | No stall |   0.157|
| YPR165W   | Stall    |   0.156|
| YGR082W   | No stall |   0.155|
| YLR354C   | No stall |   0.154|
| YPR154W   | No stall |   0.153|
| YLR344W   | No stall |   0.152|
| YDR304C   | No stall |   0.152|
| YPL213W   | No stall |   0.151|
| YPR149W   | No stall |   0.150|
| YOR123C   | Stall    |   0.149|
| YLR314C   | No stall |   0.149|
| YGL181W   | No stall |   0.149|
| YMR197C   | No stall |   0.149|
| YHL033C   | No stall |   0.148|
| YDR280W   | No stall |   0.147|
| YDR139C   | No stall |   0.147|
| YNL229C   | No stall |   0.145|
| YDR051C   | No stall |   0.145|
| YBR189W   | No stall |   0.143|
| YNR026C   | No stall |   0.143|
| YGR049W   | No stall |   0.143|
| YGR240C   | No stall |   0.142|
| YDR447C   | Stall    |   0.142|
| YNL268W   | No stall |   0.142|
| YPL106C   | No stall |   0.140|
| YLL018C   | No stall |   0.138|
| YHR141C   | No stall |   0.138|
| YGR189C   | No stall |   0.137|
| YPR102C   | Stall    |   0.136|
| YIL105C   | No stall |   0.136|
| YOR276W   | No stall |   0.135|
| YGL147C   | No stall |   0.134|
| YKL045W   | No stall |   0.134|
| YGL054C   | No stall |   0.133|
| YJL024C   | No stall |   0.132|
| YDL003W   | No stall |   0.132|
| YPR073C   | No stall |   0.131|
| YBR127C   | No stall |   0.131|
| YEL026W   | No stall |   0.130|
| YER099C   | No stall |   0.130|
| YAR002W   | No stall |   0.129|
| YBR166C   | No stall |   0.128|
| YLR250W   | No stall |   0.127|
| YJR109C   | No stall |   0.127|
| YDL092W   | Stall    |   0.127|
| YER048W-A | No stall |   0.127|
| YAL005C   | No stall |   0.126|
| YNL220W   | No stall |   0.126|
| YKR084C   | Stall    |   0.126|
| YGL040C   | No stall |   0.126|
| YMR188C   | No stall |   0.125|
| YPL265W   | No stall |   0.124|
| YLR333C   | Stall    |   0.124|
| YBR282W   | No stall |   0.124|
| YFL045C   | No stall |   0.124|
| YMR222C   | No stall |   0.123|
| YFR031C-A | No stall |   0.122|
| YOL041C   | Stall    |   0.122|
| YNL189W   | No stall |   0.122|
| YJL151C   | No stall |   0.121|
| YJL104W   | No stall |   0.120|
| YNL081C   | No stall |   0.120|
| YKR095W-A | No stall |   0.119|
| YOR230W   | No stall |   0.118|
| YNL265C   | Stall    |   0.118|
| YGR092W   | Stall    |   0.117|
| YEL027W   | No stall |   0.117|
| YHR121W   | No stall |   0.117|
| YDL191W   | Stall    |   0.117|
| YDL136W   | Stall    |   0.117|
| YFL038C   | No stall |   0.116|
| YBR191W   | No stall |   0.116|
| YDR382W   | No stall |   0.116|
| YOR202W   | No stall |   0.116|
| YMR242C   | No stall |   0.115|
| YHR005C-A | No stall |   0.112|
| YDR519W   | No stall |   0.112|
| YDR168W   | No stall |   0.111|
| YGL253W   | No stall |   0.111|
| YOR285W   | No stall |   0.110|
| YKL053C-A | No stall |   0.110|
| YNR052C   | No stall |   0.109|
| YDL046W   | No stall |   0.108|
| YMR161W   | No stall |   0.107|
| YKL180W   | No stall |   0.107|
| YBL001C   | No stall |   0.106|
| YDL097C   | No stall |   0.106|
| YLR292C   | No stall |   0.105|
| YEL054C   | Stall    |   0.105|
| YGR173W   | No stall |   0.105|
| YHR191C   | No stall |   0.104|
| YCL033C   | No stall |   0.104|
| YOL012C   | No stall |   0.103|
| YGL028C   | No stall |   0.103|
| YLR075W   | No stall |   0.102|
| YGR155W   | No stall |   0.102|
| YNL098C   | No stall |   0.102|
| YJL125C   | No stall |   0.101|
| YPL252C   | No stall |   0.101|
| YOR224C   | No stall |   0.100|
| YJL096W   | Stall    |   0.099|
| YPR069C   | No stall |   0.099|
| YJL184W   | No stall |   0.099|
| YOR344C   | Stall    |   0.098|
| YDL185W   | No stall |   0.098|
| YJL052W   | No stall |   0.097|
| YKL152C   | No stall |   0.097|
| YDL055C   | No stall |   0.096|
| YBR162W-A | Stall    |   0.096|
| YBR038W   | No stall |   0.096|
| YIL021W   | No stall |   0.095|
| YDR345C   | No stall |   0.095|
| YJR094W-A | No stall |   0.095|
| YCL055W   | No stall |   0.094|
| YGL032C   | No stall |   0.094|
| YAR007C   | No stall |   0.093|
| YJL080C   | No stall |   0.093|
| YNL069C   | No stall |   0.093|
| YIR022W   | No stall |   0.092|
| YER043C   | No stall |   0.092|
| YDR100W   | No stall |   0.091|
| YIL109C   | Stall    |   0.091|
| YNL016W   | No stall |   0.090|
| YML026C   | No stall |   0.090|
| YIL117C   | No stall |   0.090|
| YDR246W   | No stall |   0.089|
| YBR022W   | No stall |   0.089|
| YKL082C   | Stall    |   0.089|
| YDL064W   | No stall |   0.089|
| YBR214W   | No stall |   0.089|
| YOR122C   | No stall |   0.088|
| YOR157C   | No stall |   0.088|
| YDL103C   | No stall |   0.088|
| YPR032W   | No stall |   0.086|
| YML024W   | Stall    |   0.086|
| YOL094C   | No stall |   0.086|
| YER027C   | No stall |   0.086|
| YHR136C   | No stall |   0.085|
| YBL046W   | Stall    |   0.085|
| YJR007W   | No stall |   0.085|
| YDL010W   | No stall |   0.085|
| YKL016C   | No stall |   0.084|
| YKL145W   | No stall |   0.084|
| YDR511W   | No stall |   0.083|
| YNL244C   | No stall |   0.081|
| YMR079W   | No stall |   0.081|
| YBL002W   | Stall    |   0.079|
| YDR465C   | No stall |   0.078|
| YGR284C   | No stall |   0.076|
| YBL026W   | No stall |   0.075|
| YNL263C   | No stall |   0.075|
| YAL049C   | No stall |   0.075|
| YOL030W   | Stall    |   0.074|
| YMR225C   | No stall |   0.074|
| YLR259C   | No stall |   0.072|
| YGR135W   | No stall |   0.071|
| YDL133C-A | Stall    |   0.071|
| YKR059W   | Stall    |   0.071|
| YOR096W   | No stall |   0.071|
| YMR303C   | No stall |   0.071|
| YNL310C   | No stall |   0.071|
| YDL115C   | No stall |   0.070|
| YJL138C   | Stall    |   0.070|
| YBR048W   | Stall    |   0.069|
| YBR149W   | No stall |   0.069|
| YHR122W   | No stall |   0.067|
| YOR163W   | No stall |   0.067|
| YOR280C   | No stall |   0.067|
| YMR158W   | No stall |   0.067|
| YOR286W   | No stall |   0.067|
| YPL090C   | Stall    |   0.067|
| YLL024C   | No stall |   0.066|
| YBL099W   | No stall |   0.066|
| YKL060C   | No stall |   0.066|
| YDL053C   | No stall |   0.066|
| YPR080W   | No stall |   0.066|
| YJR044C   | No stall |   0.066|
| YOR182C   | Stall    |   0.065|
| YBR118W   | No stall |   0.065|
| YKL204W   | Stall    |   0.064|
| YCR031C   | No stall |   0.064|
| YMR083W   | No stall |   0.063|
| YGL123W   | No stall |   0.063|
| YOR141C   | Stall    |   0.063|
| YIL074C   | No stall |   0.062|
| YPL228W   | No stall |   0.062|
| YGL103W   | No stall |   0.062|
| YLR192C   | Stall    |   0.062|
| YNR050C   | No stall |   0.062|
| YDR450W   | No stall |   0.062|
| YOR095C   | No stall |   0.061|
| YML123C   | No stall |   0.061|
| YDL081C   | No stall |   0.061|
| YMR123W   | No stall |   0.060|
| YJL011C   | No stall |   0.060|
| YJL121C   | No stall |   0.060|
| YMR214W   | No stall |   0.060|
| YML098W   | No stall |   0.059|
| YNL071W   | No stall |   0.058|
| YOL127W   | No stall |   0.058|
| YBR039W   | No stall |   0.058|
| YPL079W   | No stall |   0.058|
| YEL071W   | No stall |   0.058|
| YMR092C   | No stall |   0.056|
| YHL001W   | Stall    |   0.056|
| YNL104C   | No stall |   0.055|
| YKL167C   | No stall |   0.055|
| YJR097W   | No stall |   0.054|
| YJL145W   | No stall |   0.054|
| YNL178W   | No stall |   0.054|
| YDR533C   | No stall |   0.053|
| YMR146C   | No stall |   0.053|
| YDR322W   | No stall |   0.053|
| YMR200W   | No stall |   0.053|
| YML080W   | No stall |   0.051|
| YER067W   | No stall |   0.050|
| YKL040C   | No stall |   0.050|
| YJL069C   | No stall |   0.049|
| YDR529C   | No stall |   0.049|
| YNL153C   | No stall |   0.048|
| YFL018C   | No stall |   0.048|
| YNL259C   | No stall |   0.047|
| YKR065C   | No stall |   0.047|
| YER026C   | No stall |   0.047|
| YPL010W   | No stall |   0.047|
| YHR010W   | Stall    |   0.046|
| YOR128C   | No stall |   0.046|
| YPL048W   | No stall |   0.046|
| YGR148C   | Stall    |   0.045|
| YFL037W   | No stall |   0.045|
| YDL130W   | No stall |   0.044|
| YDL217C   | No stall |   0.044|
| YDR233C   | No stall |   0.044|
| YEL060C   | Stall    |   0.043|
| YJR077C   | No stall |   0.043|
| YJR009C   | No stall |   0.043|
| YPL211W   | No stall |   0.043|
| YDR120C   | Stall    |   0.043|
| YOR056C   | Stall    |   0.043|
| YML022W   | No stall |   0.043|
| YOR204W   | No stall |   0.042|
| YLR136C   | No stall |   0.042|
| YGR214W   | No stall |   0.042|
| YBR160W   | No stall |   0.041|
| YPL081W   | No stall |   0.041|
| YJL063C   | No stall |   0.041|
| YOR340C   | No stall |   0.041|
| YNL075W   | No stall |   0.041|
| YOR232W   | No stall |   0.041|
| YPL086C   | No stall |   0.041|
| YBR196C   | No stall |   0.040|
| YGL056C   | No stall |   0.038|
| YDL124W   | No stall |   0.037|
| YPL139C   | No stall |   0.037|
| YFR051C   | Stall    |   0.037|
| YLR229C   | No stall |   0.036|
| YDR361C   | Stall    |   0.035|
| YIL116W   | No stall |   0.035|
| YDR063W   | No stall |   0.035|
| YLR186W   | No stall |   0.034|
| YDR454C   | No stall |   0.034|
| YPR103W   | No stall |   0.034|
| YPR016C   | No stall |   0.033|
| YBR120C   | No stall |   0.032|
| YJL146W   | No stall |   0.032|
| YPR040W   | No stall |   0.032|
| YKL007W   | No stall |   0.031|
| YBR054W   | Stall    |   0.031|
| YLL041C   | No stall |   0.030|
| YKR071C   | Stall    |   0.030|
| YDL219W   | No stall |   0.029|
| YHR107C   | No stall |   0.029|
| YHR037W   | No stall |   0.028|
| YBR121C   | Stall    |   0.028|
| YOL096C   | No stall |   0.027|
| YPR188C   | No stall |   0.027|
| YOR272W   | Stall    |   0.026|
| YOR312C   | No stall |   0.026|
| YGL145W   | No stall |   0.026|
| YGL153W   | No stall |   0.025|
| YOR046C   | No stall |   0.023|
| YIL133C   | No stall |   0.022|
| YDL165W   | No stall |   0.022|
| YOL056W   | No stall |   0.022|
| YJL173C   | No stall |   0.020|
| YHR111W   | No stall |   0.019|
| YBR101C   | No stall |   0.019|
| YDR471W   | Stall    |   0.019|
| YMR011W   | No stall |   0.019|
| YOL123W   | No stall |   0.018|
| YOL120C   | No stall |   0.018|
| YOL059W   | No stall |   0.018|
| YOL097C   | No stall |   0.018|
| YBR181C   | Stall    |   0.017|
| YER092W   | No stall |   0.017|
| YNR009W   | No stall |   0.016|
| YIL018W   | No stall |   0.015|
| YJR063W   | No stall |   0.015|
| YER133W   | No stall |   0.015|
| YBR025C   | No stall |   0.015|
| YKL214C   | No stall |   0.015|
| YKL193C   | No stall |   0.015|
| YLR048W   | No stall |   0.014|
| YMR183C   | Stall    |   0.014|
| YHR104W   | No stall |   0.014|
| YJR086W   | No stall |   0.014|
| YAL036C   | No stall |   0.014|
| YBR155W   | No stall |   0.013|
| YGR175C   | No stall |   0.013|
| YPR094W   | No stall |   0.013|
| YHR215W   | No stall |   0.013|
| YDR260C   | No stall |   0.012|
| YOR215C   | No stall |   0.012|
| YJL217W   | No stall |   0.011|
| YPL063W   | Stall    |   0.011|
| YKL196C   | No stall |   0.011|
| YPR051W   | No stall |   0.010|
| YGR074W   | No stall |   0.009|
| YMR049C   | No stall |   0.009|
| YHR076W   | No stall |   0.008|
| YKL170W   | No stall |   0.008|
| YMR208W   | No stall |   0.008|
| YGL148W   | No stall |   0.008|
| YFL016C   | No stall |   0.008|
| YPL188W   | No stall |   0.007|
| YKL046C   | No stall |   0.006|
| YKL006W   | Stall    |   0.005|
| YIL093C   | No stall |   0.002|
| YBL090W   | No stall |   0.002|
| YPR133W-A | No stall |   0.002|
| YAR071W   | No stall |   0.001|
| YML120C   | No stall |   0.000|
| YOR226C   | No stall |   0.000|
| YER141W   | No stall |  -0.001|
| YBR234C   | No stall |  -0.002|
| YGR034W   | No stall |  -0.002|
| YML086C   | Stall    |  -0.002|
| YDR175C   | Stall    |  -0.002|
| YDR087C   | No stall |  -0.002|
| YDR516C   | No stall |  -0.002|
| YGL116W   | No stall |  -0.003|
| YML030W   | No stall |  -0.004|
| YGR086C   | No stall |  -0.004|
| YFR052W   | No stall |  -0.004|
| YER095W   | No stall |  -0.005|
| YGR055W   | No stall |  -0.005|
| YLR008C   | Stall    |  -0.005|
| YDR116C   | No stall |  -0.006|
| YPL243W   | Stall    |  -0.006|
| YCR012W   | No stall |  -0.006|
| YMR142C   | Stall    |  -0.006|
| YDR298C   | No stall |  -0.007|
| YPL028W   | No stall |  -0.008|
| YGR072W   | Stall    |  -0.008|
| YMR099C   | No stall |  -0.008|
| YJL020C   | Stall    |  -0.009|
| YER165W   | No stall |  -0.009|
| YHL011C   | No stall |  -0.009|
| YOR310C   | Stall    |  -0.009|
| YKL127W   | No stall |  -0.010|
| YBL087C   | Stall    |  -0.010|
| YPL196W   | No stall |  -0.010|
| YOR117W   | No stall |  -0.011|
| YOR323C   | No stall |  -0.011|
| YKR055W   | No stall |  -0.012|
| YGR207C   | No stall |  -0.012|
| YHR013C   | No stall |  -0.012|
| YHL027W   | No stall |  -0.013|
| YKL021C   | Stall    |  -0.013|
| YKL190W   | No stall |  -0.014|
| YOR362C   | No stall |  -0.014|
| YMR022W   | No stall |  -0.014|
| YBR227C   | No stall |  -0.014|
| YPR131C   | No stall |  -0.015|
| YDR418W   | Stall    |  -0.015|
| YNR003C   | No stall |  -0.016|
| YHR049W   | No stall |  -0.017|
| YDR494W   | No stall |  -0.017|
| YHR208W   | No stall |  -0.018|
| YIL142W   | No stall |  -0.019|
| YHR128W   | No stall |  -0.020|
| YML048W   | No stall |  -0.020|
| YOR184W   | No stall |  -0.022|
| YKL099C   | Stall    |  -0.022|
| YGL231C   | No stall |  -0.022|
| YCR002C   | No stall |  -0.023|
| YKL029C   | No stall |  -0.024|
| YDR486C   | No stall |  -0.025|
| YDL126C   | No stall |  -0.025|
| YMR153W   | No stall |  -0.026|
| YBR268W   | Stall    |  -0.027|
| YJL177W   | No stall |  -0.027|
| YML060W   | No stall |  -0.027|
| YHR025W   | No stall |  -0.028|
| YGR038W   | No stall |  -0.028|
| YGR081C   | Stall    |  -0.029|
| YBR091C   | No stall |  -0.029|
| YLL027W   | Stall    |  -0.030|
| YJR088C   | No stall |  -0.030|
| YPL218W   | No stall |  -0.031|
| YDR451C   | Stall    |  -0.032|
| YOR347C   | No stall |  -0.032|
| YKL025C   | Stall    |  -0.032|
| YKL141W   | No stall |  -0.032|
| YLR027C   | No stall |  -0.033|
| YLR369W   | No stall |  -0.033|
| YNL222W   | No stall |  -0.033|
| YDR025W   | Stall    |  -0.033|
| YIL016W   | Stall    |  -0.034|
| YBR221C   | No stall |  -0.034|
| YMR105C   | No stall |  -0.034|
| YGL020C   | No stall |  -0.034|
| YGR106C   | No stall |  -0.034|
| YLR347C   | No stall |  -0.034|
| YGL225W   | No stall |  -0.035|
| YML025C   | No stall |  -0.036|
| YLR199C   | No stall |  -0.038|
| YBL091C-A | No stall |  -0.038|
| YDR394W   | No stall |  -0.039|
| YGR159C   | No stall |  -0.041|
| YFR047C   | No stall |  -0.041|
| YHR207C   | No stall |  -0.041|
| YER136W   | No stall |  -0.042|
| YPL220W   | Stall    |  -0.043|
| YGL122C   | No stall |  -0.043|
| YGR253C   | No stall |  -0.043|
| YEL015W   | No stall |  -0.043|
| YGR060W   | No stall |  -0.044|
| YOL139C   | No stall |  -0.044|
| YKL088W   | No stall |  -0.045|
| YDR321W   | No stall |  -0.045|
| YHR167W   | No stall |  -0.045|
| YLR216C   | No stall |  -0.046|
| YBR066C   | Stall    |  -0.046|
| YIL004C   | No stall |  -0.046|
| YJR118C   | No stall |  -0.047|
| YCR009C   | No stall |  -0.047|
| YPR133C   | Stall    |  -0.048|
| YLR017W   | No stall |  -0.049|
| YDL215C   | No stall |  -0.050|
| YPL154C   | No stall |  -0.050|
| YIL065C   | No stall |  -0.051|
| YDR192C   | No stall |  -0.051|
| YDR121W   | No stall |  -0.053|
| YNL112W   | No stall |  -0.053|
| YDL137W   | No stall |  -0.054|
| YGL035C   | Stall    |  -0.054|
| YDR151C   | Stall    |  -0.054|
| YGL068W   | No stall |  -0.056|
| YOR004W   | Stall    |  -0.056|
| YDR214W   | No stall |  -0.056|
| YMR108W   | No stall |  -0.057|
| YGL189C   | No stall |  -0.057|
| YFR024C-A | Stall    |  -0.057|
| YER044C   | No stall |  -0.057|
| YHR187W   | No stall |  -0.057|
| YPL055C   | No stall |  -0.058|
| YPR033C   | No stall |  -0.058|
| YGL105W   | No stall |  -0.059|
| YMR152W   | No stall |  -0.059|
| YDL202W   | No stall |  -0.059|
| YOL064C   | No stall |  -0.059|
| YHR060W   | No stall |  -0.061|
| YMR205C   | No stall |  -0.061|
| YFR044C   | No stall |  -0.061|
| YPR129W   | No stall |  -0.062|
| YOR194C   | No stall |  -0.064|
| YDR083W   | Stall    |  -0.064|
| YBL038W   | No stall |  -0.065|
| YHR019C   | No stall |  -0.065|
| YER103W   | No stall |  -0.066|
| YDR358W   | No stall |  -0.066|
| YBR052C   | No stall |  -0.066|
| YIL083C   | No stall |  -0.066|
| YML004C   | No stall |  -0.067|
| YKL120W   | No stall |  -0.067|
| YNL207W   | Stall    |  -0.067|
| YER146W   | No stall |  -0.067|
| YDR405W   | No stall |  -0.068|
| YKR093W   | No stall |  -0.070|
| YJL036W   | No stall |  -0.070|
| YIL076W   | No stall |  -0.071|
| YPR145W   | No stall |  -0.072|
| YJL174W   | No stall |  -0.072|
| YDR339C   | No stall |  -0.072|
| YLR193C   | No stall |  -0.072|
| YDL198C   | No stall |  -0.073|
| YDL226C   | No stall |  -0.073|
| YIL022W   | Stall    |  -0.073|
| YBR037C   | No stall |  -0.074|
| YIL157C   | No stall |  -0.074|
| YML100W   | No stall |  -0.075|
| YDR054C   | No stall |  -0.075|
| YNR035C   | No stall |  -0.076|
| YDR153C   | Stall    |  -0.076|
| YOR176W   | No stall |  -0.078|
| YCL059C   | Stall    |  -0.078|
| YNR018W   | No stall |  -0.079|
| YJR133W   | No stall |  -0.079|
| YDR378C   | No stall |  -0.079|
| YML102W   | No stall |  -0.080|
| YDL201W   | Stall    |  -0.080|
| YER131W   | No stall |  -0.080|
| YLR304C   | No stall |  -0.080|
| YML057W   | No stall |  -0.080|
| YBR263W   | No stall |  -0.081|
| YKR066C   | No stall |  -0.081|
| YPL169C   | No stall |  -0.081|
| YHR064C   | No stall |  -0.082|
| YLR060W   | No stall |  -0.082|
| YBR034C   | No stall |  -0.083|
| YGL161C   | No stall |  -0.083|
| YBR170C   | No stall |  -0.083|
| YGL135W   | Stall    |  -0.084|
| YER091C   | No stall |  -0.084|
| YMR120C   | No stall |  -0.085|
| YMR116C   | No stall |  -0.085|
| YCR072C   | No stall |  -0.085|
| YDL227C   | No stall |  -0.086|
| YMR281W   | No stall |  -0.086|
| YOR181W   | Stall    |  -0.086|
| YMR009W   | No stall |  -0.087|
| YIL130W   | No stall |  -0.087|
| YPR191W   | No stall |  -0.087|
| YER110C   | No stall |  -0.087|
| YNL138W   | Stall    |  -0.087|
| YMR043W   | No stall |  -0.088|
| YER004W   | No stall |  -0.088|
| YGL011C   | No stall |  -0.088|
| YNL248C   | No stall |  -0.088|
| YGR211W   | No stall |  -0.088|
| YGR030C   | No stall |  -0.088|
| YJR093C   | Stall    |  -0.089|
| YFL028C   | No stall |  -0.089|
| YMR195W   | No stall |  -0.089|
| YLL010C   | No stall |  -0.090|
| YGL210W   | No stall |  -0.090|
| YJL010C   | Stall    |  -0.090|
| YGL048C   | No stall |  -0.090|
| YBL058W   | No stall |  -0.091|
| YGL098W   | No stall |  -0.091|
| YHR132W-A | No stall |  -0.091|
| YIL010W   | No stall |  -0.091|
| YHR070W   | No stall |  -0.092|
| YJL041W   | No stall |  -0.092|
| YDR183W   | No stall |  -0.092|
| YGR282C   | No stall |  -0.093|
| YOR164C   | No stall |  -0.094|
| YOL027C   | No stall |  -0.094|
| YOL018C   | No stall |  -0.095|
| YER082C   | No stall |  -0.095|
| YJR072C   | Stall    |  -0.097|
| YNR049C   | No stall |  -0.097|
| YJL111W   | No stall |  -0.099|
| YOR251C   | No stall |  -0.100|
| YER012W   | No stall |  -0.100|
| YKL143W   | Stall    |  -0.101|
| YKL157W   | No stall |  -0.101|
| YEL003W   | No stall |  -0.102|
| YDR188W   | Stall    |  -0.102|
| YIL125W   | No stall |  -0.102|
| YIL063C   | No stall |  -0.102|
| YKL013C   | No stall |  -0.103|
| YER036C   | Stall    |  -0.103|
| YGL220W   | No stall |  -0.103|
| YKL148C   | No stall |  -0.103|
| YER056C   | No stall |  -0.104|
| YLR359W   | No stall |  -0.104|
| YFR037C   | No stall |  -0.104|
| YJL078C   | No stall |  -0.104|
| YIL041W   | No stall |  -0.105|
| YPL203W   | No stall |  -0.105|
| YDR429C   | No stall |  -0.105|
| YCL040W   | No stall |  -0.106|
| YBL040C   | No stall |  -0.106|
| YBR165W   | No stall |  -0.107|
| YPR187W   | No stall |  -0.108|
| YOR145C   | No stall |  -0.108|
| YPR124W   | No stall |  -0.110|
| YJR148W   | No stall |  -0.110|
| YGR203W   | No stall |  -0.110|
| YFR010W   | Stall    |  -0.110|
| YOL007C   | No stall |  -0.110|
| YGL157W   | No stall |  -0.110|
| YJL097W   | No stall |  -0.111|
| YPL031C   | No stall |  -0.111|
| YPL173W   | No stall |  -0.111|
| YPR074C   | No stall |  -0.111|
| YPL249C   | No stall |  -0.111|
| YNL241C   | No stall |  -0.111|
| YLR120C   | No stall |  -0.112|
| YMR309C   | No stall |  -0.112|
| YPR180W   | No stall |  -0.112|
| YPR107C   | No stall |  -0.113|
| YBL032W   | No stall |  -0.114|
| YGL091C   | No stall |  -0.114|
| YPR112C   | Stall    |  -0.114|
| YGL213C   | No stall |  -0.114|
| YIL038C   | No stall |  -0.114|
| YPL099C   | No stall |  -0.115|
| YOR185C   | No stall |  -0.115|
| YBL047C   | No stall |  -0.116|
| YML017W   | No stall |  -0.116|
| YPL190C   | Stall    |  -0.117|
| YOL147C   | No stall |  -0.118|
| YBR111W-A | No stall |  -0.118|
| YBR273C   | No stall |  -0.119|
| YKL056C   | No stall |  -0.119|
| YLR286C   | No stall |  -0.119|
| YJR099W   | No stall |  -0.119|
| YLR002C   | Stall    |  -0.120|
| YIR038C   | No stall |  -0.120|
| YLR195C   | No stall |  -0.121|
| YOR061W   | No stall |  -0.121|
| YDL164C   | No stall |  -0.121|
| YKL211C   | No stall |  -0.121|
| YDR388W   | No stall |  -0.122|
| YDR308C   | No stall |  -0.122|
| YEL063C   | No stall |  -0.122|
| YOL052C   | No stall |  -0.122|
| YBR176W   | No stall |  -0.122|
| YKL210W   | No stall |  -0.123|
| YMR136W   | Stall    |  -0.123|
| YJL066C   | No stall |  -0.123|
| YGR172C   | No stall |  -0.124|
| YHR063C   | No stall |  -0.124|
| YDR337W   | Stall    |  -0.124|
| YNL074C   | No stall |  -0.124|
| YJL030W   | No stall |  -0.124|
| YKL138C   | Stall    |  -0.124|
| YHR200W   | No stall |  -0.125|
| YFR014C   | No stall |  -0.125|
| YPR098C   | No stall |  -0.126|
| YER145C   | No stall |  -0.127|
| YLR118C   | No stall |  -0.127|
| YMR300C   | No stall |  -0.127|
| YCL001W   | No stall |  -0.128|
| YGL234W   | No stall |  -0.128|
| YGL209W   | No stall |  -0.128|
| YHR001W   | No stall |  -0.128|
| YDL193W   | No stall |  -0.129|
| YHR029C   | No stall |  -0.129|
| YOR065W   | No stall |  -0.129|
| YBL027W   | No stall |  -0.130|
| YPR181C   | No stall |  -0.131|
| YGR208W   | No stall |  -0.131|
| YLR268W   | No stall |  -0.131|
| YML014W   | No stall |  -0.132|
| YDR068W   | No stall |  -0.132|
| YOR069W   | No stall |  -0.132|
| YLR356W   | No stall |  -0.132|
| YKL144C   | No stall |  -0.133|
| YHR158C   | No stall |  -0.133|
| YDL100C   | No stall |  -0.133|
| YDR091C   | No stall |  -0.134|
| YHR027C   | No stall |  -0.134|
| YKL024C   | No stall |  -0.134|
| YPL158C   | Stall    |  -0.136|
| YDR341C   | No stall |  -0.136|
| YGR020C   | No stall |  -0.137|
| YNL141W   | No stall |  -0.140|
| YDL015C   | Stall    |  -0.140|
| YLR433C   | No stall |  -0.140|
| YLR231C   | No stall |  -0.141|
| YDR084C   | No stall |  -0.141|
| YJL178C   | No stall |  -0.141|
| YLR285W   | No stall |  -0.141|
| YCL030C   | No stall |  -0.142|
| YPL078C   | No stall |  -0.143|
| YDR477W   | No stall |  -0.144|
| YDR033W   | Stall    |  -0.144|
| YLR093C   | No stall |  -0.144|
| YGR095C   | No stall |  -0.145|
| YDL084W   | No stall |  -0.146|
| YER094C   | No stall |  -0.147|
| YPL146C   | Stall    |  -0.147|
| YBR143C   | No stall |  -0.147|
| YPR100W   | No stall |  -0.147|
| YDR075W   | No stall |  -0.147|
| YPR198W   | No stall |  -0.148|
| YBR077C   | No stall |  -0.148|
| YJL060W   | No stall |  -0.149|
| YDL229W   | No stall |  -0.151|
| YPL181W   | Stall    |  -0.151|
| YPR156C   | Stall    |  -0.152|
| YBR261C   | No stall |  -0.155|
| YPL237W   | Stall    |  -0.156|
| YLR148W   | No stall |  -0.156|
| YML085C   | No stall |  -0.156|
| YBR084C-A | No stall |  -0.157|
| YDL131W   | No stall |  -0.157|
| YFR003C   | Stall    |  -0.157|
| YMR305C   | No stall |  -0.158|
| YJL191W   | Stall    |  -0.159|
| YKR048C   | No stall |  -0.159|
| YGL047W   | No stall |  -0.160|
| YMR217W   | No stall |  -0.160|
| YGR199W   | No stall |  -0.160|
| YDR487C   | No stall |  -0.160|
| YMR110C   | No stall |  -0.160|
| YMR015C   | No stall |  -0.161|
| YGL050W   | Stall    |  -0.161|
| YPR110C   | No stall |  -0.163|
| YOR359W   | No stall |  -0.163|
| YNL027W   | Stall    |  -0.165|
| YOL119C   | No stall |  -0.165|
| YHR113W   | No stall |  -0.165|
| YOL010W   | No stall |  -0.165|
| YGL221C   | No stall |  -0.165|
| YGR158C   | No stall |  -0.166|
| YBL064C   | No stall |  -0.167|
| YIR008C   | No stall |  -0.167|
| YML062C   | Stall    |  -0.167|
| YBL093C   | No stall |  -0.167|
| YGL009C   | No stall |  -0.169|
| YNL209W   | No stall |  -0.169|
| YDR496C   | No stall |  -0.169|
| YOL016C   | No stall |  -0.170|
| YBR171W   | No stall |  -0.171|
| YPL204W   | No stall |  -0.172|
| YNR041C   | No stall |  -0.173|
| YOL026C   | No stall |  -0.173|
| YLR441C   | Stall    |  -0.174|
| YNL116W   | No stall |  -0.174|
| YER049W   | No stall |  -0.175|
| YNL212W   | No stall |  -0.175|
| YHR143W   | No stall |  -0.175|
| YHR058C   | No stall |  -0.175|
| YDR097C   | No stall |  -0.176|
| YOR361C   | Stall    |  -0.176|
| YHR057C   | No stall |  -0.177|
| YNR032C-A | No stall |  -0.177|
| YDL178W   | No stall |  -0.177|
| YLR209C   | No stall |  -0.178|
| YCL054W   | Stall    |  -0.179|
| YJL014W   | No stall |  -0.179|
| YDR472W   | No stall |  -0.179|
| YER023W   | No stall |  -0.180|
| YGL025C   | Stall    |  -0.180|
| YAL016W   | No stall |  -0.180|
| YKL003C   | No stall |  -0.180|
| YKR092C   | No stall |  -0.180|
| YBR198C   | No stall |  -0.180|
| YHL002W   | No stall |  -0.181|
| YKL130C   | No stall |  -0.181|
| YBR080C   | No stall |  -0.182|
| YOL070C   | Stall    |  -0.182|
| YML013W   | No stall |  -0.182|
| YLR220W   | No stall |  -0.183|
| YGL001C   | No stall |  -0.183|
| YER107C   | No stall |  -0.183|
| YIL124W   | No stall |  -0.183|
| YLR178C   | No stall |  -0.184|
| YDR036C   | No stall |  -0.184|
| YOR253W   | No stall |  -0.184|
| YMR290C   | No stall |  -0.185|
| YHR127W   | No stall |  -0.185|
| YGR124W   | No stall |  -0.186|
| YHL004W   | No stall |  -0.187|
| YLR044C   | No stall |  -0.187|
| YOL135C   | No stall |  -0.188|
| YER037W   | Stall    |  -0.189|
| YOR091W   | Stall    |  -0.190|
| YIR004W   | No stall |  -0.190|
| YOL129W   | No stall |  -0.191|
| YOR168W   | Stall    |  -0.191|
| YCR071C   | No stall |  -0.192|
| YLL001W   | No stall |  -0.192|
| YGL186C   | Stall    |  -0.192|
| YJL003W   | No stall |  -0.193|
| YDR258C   | No stall |  -0.194|
| YDR281C   | No stall |  -0.195|
| YNL316C   | No stall |  -0.195|
| YKL185W   | Stall    |  -0.196|
| YOL061W   | No stall |  -0.198|
| YPR019W   | No stall |  -0.200|
| YOR288C   | No stall |  -0.200|
| YDL047W   | No stall |  -0.202|
| YBR159W   | No stall |  -0.203|
| YGL200C   | No stall |  -0.203|
| YGR116W   | Stall    |  -0.203|
| YDL040C   | No stall |  -0.203|
| YNL252C   | No stall |  -0.204|
| YOR085W   | No stall |  -0.205|
| YDL147W   | No stall |  -0.205|
| YPR176C   | No stall |  -0.206|
| YNR001C   | No stall |  -0.206|
| YDL108W   | No stall |  -0.206|
| YBR024W   | Stall    |  -0.207|
| YKL186C   | No stall |  -0.208|
| YJL122W   | Stall    |  -0.208|
| YDL182W   | No stall |  -0.209|
| YEL012W   | No stall |  -0.209|
| YIL034C   | No stall |  -0.210|
| YOL142W   | No stall |  -0.210|
| YNL177C   | No stall |  -0.211|
| YER068W   | Stall    |  -0.211|
| YOL126C   | No stall |  -0.212|
| YJL087C   | No stall |  -0.212|
| YDR035W   | Stall    |  -0.212|
| YGR244C   | No stall |  -0.213|
| YCR035C   | No stall |  -0.213|
| YML021C   | No stall |  -0.214|
| YDR330W   | No stall |  -0.214|
| YPR028W   | No stall |  -0.215|
| YKR083C   | No stall |  -0.216|
| YEL051W   | No stall |  -0.216|
| YMR278W   | No stall |  -0.216|
| YMR215W   | No stall |  -0.216|
| YER059W   | No stall |  -0.216|
| YGR205W   | No stall |  -0.217|
| YBR177C   | No stall |  -0.218|
| YJR117W   | No stall |  -0.218|
| YNL044W   | No stall |  -0.218|
| YNL289W   | No stall |  -0.218|
| YNL132W   | Stall    |  -0.218|
| YOL022C   | No stall |  -0.219|
| YJR065C   | No stall |  -0.220|
| YOR042W   | Stall    |  -0.220|
| YLR432W   | No stall |  -0.221|
| YGR266W   | No stall |  -0.221|
| YJL008C   | No stall |  -0.221|
| YOL077C   | No stall |  -0.222|
| YKL027W   | Stall    |  -0.222|
| YFR034C   | No stall |  -0.225|
| YLL038C   | No stall |  -0.225|
| YHR123W   | No stall |  -0.225|
| YDR346C   | No stall |  -0.225|
| YOR171C   | No stall |  -0.226|
| YKR026C   | No stall |  -0.226|
| YHR016C   | Stall    |  -0.226|
| YAR002C-A | No stall |  -0.227|
| YDR408C   | No stall |  -0.227|
| YDL006W   | No stall |  -0.227|
| YOR245C   | No stall |  -0.227|
| YPL262W   | No stall |  -0.228|
| YBL023C   | Stall    |  -0.229|
| YDR189W   | No stall |  -0.229|
| YJR032W   | No stall |  -0.229|
| YBR112C   | No stall |  -0.230|
| YDL072C   | No stall |  -0.231|
| YMR067C   | No stall |  -0.232|
| YPL160W   | No stall |  -0.232|
| YCR084C   | No stall |  -0.234|
| YNL121C   | No stall |  -0.234|
| YMR241W   | No stall |  -0.234|
| YOL124C   | No stall |  -0.234|
| YCL017C   | No stall |  -0.235|
| YGR218W   | No stall |  -0.235|
| YPR060C   | No stall |  -0.236|
| YER070W   | No stall |  -0.236|
| YKR068C   | No stall |  -0.236|
| YPL274W   | No stall |  -0.237|
| YOR016C   | No stall |  -0.238|
| YOR369C   | No stall |  -0.239|
| YNL251C   | Stall    |  -0.239|
| YCR026C   | No stall |  -0.239|
| YGR262C   | No stall |  -0.239|
| YDR399W   | No stall |  -0.240|
| YJR073C   | No stall |  -0.240|
| YKR008W   | No stall |  -0.240|
| YFR009W   | No stall |  -0.241|
| YKL212W   | No stall |  -0.242|
| YLR420W   | No stall |  -0.242|
| YJR074W   | No stall |  -0.243|
| YGL084C   | No stall |  -0.244|
| YDL066W   | No stall |  -0.244|
| YIL061C   | No stall |  -0.245|
| YDL095W   | No stall |  -0.245|
| YBL039C   | No stall |  -0.246|
| YNL232W   | No stall |  -0.247|
| YOR008C   | No stall |  -0.248|
| YKL081W   | No stall |  -0.248|
| YJL171C   | No stall |  -0.248|
| YMR033W   | No stall |  -0.248|
| YNL308C   | Stall    |  -0.248|
| YDR404C   | No stall |  -0.249|
| YML055W   | No stall |  -0.250|
| YDR021W   | No stall |  -0.250|
| YPR161C   | No stall |  -0.251|
| YER119C   | No stall |  -0.251|
| YKL069W   | No stall |  -0.251|
| YDL232W   | No stall |  -0.251|
| YPL232W   | Stall    |  -0.252|
| YIL078W   | No stall |  -0.252|
| YGL171W   | Stall    |  -0.253|
| YDR001C   | No stall |  -0.253|
| YEL058W   | No stall |  -0.254|
| YDL045C   | No stall |  -0.254|
| YOL149W   | No stall |  -0.256|
| YOR316C   | No stall |  -0.256|
| YAL030W   | No stall |  -0.256|
| YDR284C   | No stall |  -0.257|
| YPR034W   | No stall |  -0.257|
| YAL046C   | No stall |  -0.257|
| YDR270W   | No stall |  -0.257|
| YKR074W   | No stall |  -0.258|
| YGR054W   | Stall    |  -0.258|
| YIL040W   | No stall |  -0.258|
| YEL001C   | No stall |  -0.258|
| YBR257W   | Stall    |  -0.259|
| YMR314W   | No stall |  -0.259|
| YOR212W   | No stall |  -0.259|
| YCL043C   | No stall |  -0.259|
| YLR332W   | No stall |  -0.259|
| YIL135C   | No stall |  -0.260|
| YKL119C   | Stall    |  -0.262|
| YPR018W   | Stall    |  -0.262|
| YJL134W   | No stall |  -0.263|
| YDR390C   | No stall |  -0.263|
| YJL124C   | No stall |  -0.264|
| YGR162W   | Stall    |  -0.264|
| YHR216W   | No stall |  -0.265|
| YLL014W   | No stall |  -0.265|
| YNL309W   | No stall |  -0.265|
| YMR004W   | No stall |  -0.265|
| YHR147C   | No stall |  -0.266|
| YBR122C   | Stall    |  -0.266|
| YGR174C   | No stall |  -0.266|
| YIL027C   | No stall |  -0.267|
| YPL145C   | No stall |  -0.267|
| YLR312W-A | No stall |  -0.267|
| YIL145C   | No stall |  -0.268|
| YEL047C   | No stall |  -0.268|
| YMR272C   | No stall |  -0.269|
| YPL065W   | No stall |  -0.269|
| YML016C   | No stall |  -0.270|
| YBR102C   | No stall |  -0.270|
| YLR183C   | No stall |  -0.270|
| YGL019W   | No stall |  -0.270|
| YKR046C   | No stall |  -0.271|
| YKL084W   | No stall |  -0.272|
| YML130C   | No stall |  -0.272|
| YER122C   | No stall |  -0.273|
| YKL139W   | No stall |  -0.273|
| YOR126C   | No stall |  -0.273|
| YMR062C   | No stall |  -0.273|
| YBR279W   | Stall    |  -0.274|
| YHR203C   | No stall |  -0.274|
| YER176W   | Stall    |  -0.275|
| YDR300C   | No stall |  -0.275|
| YPL266W   | No stall |  -0.276|
| YDR016C   | No stall |  -0.276|
| YEL017W   | No stall |  -0.276|
| YER118C   | No stall |  -0.277|
| YML052W   | No stall |  -0.278|
| YDL076C   | No stall |  -0.278|
| YBR104W   | No stall |  -0.279|
| YOR273C   | Stall    |  -0.279|
| YIL137C   | No stall |  -0.279|
| YDR373W   | No stall |  -0.279|
| YNL185C   | No stall |  -0.280|
| YBL054W   | No stall |  -0.280|
| YPR159W   | No stall |  -0.282|
| YGR185C   | Stall    |  -0.282|
| YKL104C   | No stall |  -0.282|
| YML036W   | No stall |  -0.282|
| YGR264C   | No stall |  -0.282|
| YNL288W   | No stall |  -0.282|
| YEL038W   | No stall |  -0.284|
| YGR119C   | No stall |  -0.284|
| YKR080W   | No stall |  -0.284|
| YNR037C   | No stall |  -0.285|
| YCR036W   | No stall |  -0.287|
| YCR034W   | No stall |  -0.287|
| YJR092W   | No stall |  -0.287|
| YOR179C   | No stall |  -0.287|
| YNL216W   | No stall |  -0.287|
| YKL001C   | No stall |  -0.288|
| YGR195W   | No stall |  -0.288|
| YNR053C   | Stall    |  -0.288|
| YDR195W   | No stall |  -0.288|
| YPR036W   | No stall |  -0.289|
| YLR175W   | Stall    |  -0.289|
| YOL146W   | No stall |  -0.290|
| YGR163W   | No stall |  -0.291|
| YNL192W   | No stall |  -0.291|
| YGL245W   | No stall |  -0.292|
| YDL018C   | Stall    |  -0.292|
| YMR236W   | No stall |  -0.293|
| YML010W   | Stall    |  -0.293|
| YKL126W   | No stall |  -0.293|
| YFR041C   | No stall |  -0.293|
| YIL015W   | No stall |  -0.293|
| YBR192W   | No stall |  -0.295|
| YBR185C   | No stall |  -0.295|
| YPL234C   | No stall |  -0.295|
| YLR340W   | No stall |  -0.296|
| YPL240C   | Stall    |  -0.297|
| YDR229W   | No stall |  -0.298|
| YCR020C-A | No stall |  -0.298|
| YAR027W   | No stall |  -0.299|
| YGL130W   | No stall |  -0.299|
| YGL232W   | No stall |  -0.299|
| YGL194C   | No stall |  -0.299|
| YKL018W   | No stall |  -0.299|
| YLL062C   | No stall |  -0.300|
| YPL043W   | Stall    |  -0.300|
| YMR074C   | No stall |  -0.300|
| YDR041W   | No stall |  -0.300|
| YLR390W   | Stall    |  -0.301|
| YPL051W   | No stall |  -0.302|
| YBR106W   | No stall |  -0.303|
| YPL178W   | No stall |  -0.303|
| YNR022C   | No stall |  -0.304|
| YGL127C   | No stall |  -0.304|
| YKR029C   | Stall    |  -0.304|
| YOR335C   | No stall |  -0.305|
| YLR204W   | Stall    |  -0.305|
| YLR380W   | No stall |  -0.305|
| YGR084C   | No stall |  -0.306|
| YIR011C   | No stall |  -0.306|
| YDL190C   | No stall |  -0.306|
| YMR170C   | No stall |  -0.307|
| YHR012W   | No stall |  -0.307|
| YOL004W   | No stall |  -0.308|
| YBR254C   | No stall |  -0.308|
| YOR261C   | No stall |  -0.308|
| YJL072C   | No stall |  -0.309|
| YMR024W   | No stall |  -0.309|
| YBR036C   | No stall |  -0.309|
| YLR094C   | Stall    |  -0.309|
| YOL137W   | No stall |  -0.309|
| YLR335W   | Stall    |  -0.310|
| YHR043C   | No stall |  -0.310|
| YHR142W   | No stall |  -0.310|
| YBR207W   | No stall |  -0.310|
| YJR145C   | No stall |  -0.311|
| YLR387C   | Stall    |  -0.311|
| YOR187W   | No stall |  -0.311|
| YKL002W   | No stall |  -0.311|
| YPR025C   | No stall |  -0.311|
| YDR320C   | No stall |  -0.312|
| YBR111C   | No stall |  -0.312|
| YDR320C-A | No stall |  -0.312|
| YLR181C   | No stall |  -0.312|
| YLR191W   | No stall |  -0.314|
| YML041C   | No stall |  -0.314|
| YML093W   | Stall    |  -0.315|
| YNL239W   | No stall |  -0.315|
| YNL284C   | No stall |  -0.315|
| YPR079W   | No stall |  -0.317|
| YOR264W   | No stall |  -0.317|
| YOL102C   | No stall |  -0.317|
| YER020W   | No stall |  -0.318|
| YDR237W   | No stall |  -0.318|
| YNL290W   | No stall |  -0.318|
| YER148W   | No stall |  -0.319|
| YNL084C   | Stall    |  -0.319|
| YIR003W   | Stall    |  -0.319|
| YDR201W   | Stall    |  -0.319|
| YMR014W   | Stall    |  -0.320|
| YHR201C   | No stall |  -0.320|
| YLR172C   | No stall |  -0.320|
| YGR128C   | No stall |  -0.321|
| YMR238W   | No stall |  -0.322|
| YIR026C   | No stall |  -0.322|
| YLR034C   | No stall |  -0.322|
| YDL192W   | No stall |  -0.322|
| YNL107W   | No stall |  -0.322|
| YDR057W   | No stall |  -0.322|
| YLL022C   | No stall |  -0.322|
| YML012W   | No stall |  -0.324|
| YNR017W   | No stall |  -0.324|
| YGR156W   | Stall    |  -0.326|
| YOR101W   | No stall |  -0.326|
| YMR255W   | No stall |  -0.326|
| YHR069C   | No stall |  -0.328|
| YJL196C   | No stall |  -0.328|
| YKL184W   | No stall |  -0.330|
| YER022W   | No stall |  -0.330|
| YNL317W   | No stall |  -0.331|
| YGL049C   | Stall    |  -0.331|
| YCL032W   | No stall |  -0.331|
| YEL002C   | No stall |  -0.332|
| YHL003C   | No stall |  -0.332|
| YJL172W   | No stall |  -0.332|
| YDR463W   | Stall    |  -0.332|
| YGL008C   | No stall |  -0.333|
| YKR060W   | No stall |  -0.333|
| YNL001W   | No stall |  -0.334|
| YLR429W   | No stall |  -0.334|
| YNL264C   | No stall |  -0.335|
| YLR028C   | No stall |  -0.335|
| YHR068W   | No stall |  -0.336|
| YAL007C   | No stall |  -0.336|
| YIL064W   | No stall |  -0.336|
| YDL235C   | No stall |  -0.336|
| YKR030W   | No stall |  -0.337|
| YNR036C   | No stall |  -0.337|
| YMR260C   | Stall    |  -0.338|
| YNL330C   | No stall |  -0.338|
| YNL246W   | No stall |  -0.338|
| YBR264C   | No stall |  -0.338|
| YPR072W   | Stall    |  -0.339|
| YPR097W   | Stall    |  -0.339|
| YLR323C   | No stall |  -0.340|
| YGR268C   | No stall |  -0.341|
| YMR246W   | No stall |  -0.341|
| YNL243W   | No stall |  -0.342|
| YGR142W   | No stall |  -0.342|
| YLR100W   | No stall |  -0.342|
| YNL062C   | No stall |  -0.344|
| YBR057C   | No stall |  -0.344|
| YPR086W   | No stall |  -0.344|
| YOR209C   | No stall |  -0.345|
| YHR117W   | Stall    |  -0.345|
| YCL009C   | No stall |  -0.346|
| YGR232W   | No stall |  -0.347|
| YMR301C   | No stall |  -0.348|
| YLR309C   | Stall    |  -0.350|
| YKL065C   | No stall |  -0.351|
| YJR017C   | No stall |  -0.351|
| YNL045W   | No stall |  -0.351|
| YER080W   | No stall |  -0.352|
| YJL012C   | No stall |  -0.353|
| YML035C   | No stall |  -0.353|
| YNL032W   | No stall |  -0.354|
| YCR005C   | No stall |  -0.354|
| YOL017W   | Stall    |  -0.354|
| YOR079C   | No stall |  -0.356|
| YDR395W   | No stall |  -0.357|
| YOL080C   | No stall |  -0.357|
| YLR300W   | No stall |  -0.357|
| YGL055W   | Stall    |  -0.357|
| YGL097W   | No stall |  -0.358|
| YLR134W   | No stall |  -0.359|
| YOL008W   | No stall |  -0.359|
| YGR132C   | No stall |  -0.359|
| YPL177C   | No stall |  -0.359|
| YIL036W   | No stall |  -0.361|
| YBR139W   | Stall    |  -0.361|
| YLR144C   | No stall |  -0.362|
| YNL281W   | No stall |  -0.363|
| YHR047C   | No stall |  -0.363|
| YGR277C   | No stall |  -0.363|
| YCL025C   | No stall |  -0.364|
| YEL018W   | No stall |  -0.364|
| YPR037C   | No stall |  -0.364|
| YGR090W   | No stall |  -0.364|
| YGR036C   | No stall |  -0.364|
| YEL066W   | No stall |  -0.365|
| YDR019C   | No stall |  -0.365|
| YHR149C   | No stall |  -0.367|
| YDR240C   | No stall |  -0.368|
| YNL161W   | No stall |  -0.368|
| YAL060W   | No stall |  -0.368|
| YBL103C   | Stall    |  -0.369|
| YDR453C   | No stall |  -0.370|
| YCR052W   | No stall |  -0.371|
| YDR205W   | No stall |  -0.371|
| YBR094W   | No stall |  -0.372|
| YHR100C   | No stall |  -0.372|
| YJL208C   | No stall |  -0.373|
| YOL057W   | No stall |  -0.373|
| YNL247W   | No stall |  -0.373|
| YKL110C   | No stall |  -0.374|
| YDL120W   | No stall |  -0.375|
| YLR447C   | No stall |  -0.375|
| YDR518W   | No stall |  -0.376|
| YOR143C   | No stall |  -0.376|
| YPL231W   | No stall |  -0.376|
| YCR083W   | No stall |  -0.377|
| YBL016W   | No stall |  -0.377|
| YGL099W   | Stall    |  -0.377|
| YDR047W   | No stall |  -0.377|
| YHR006W   | Stall    |  -0.378|
| YHR137W   | No stall |  -0.378|
| YCR077C   | Stall    |  -0.379|
| YJR103W   | No stall |  -0.379|
| YNR016C   | No stall |  -0.379|
| YGR267C   | No stall |  -0.379|
| YNL173C   | Stall    |  -0.380|
| YKL091C   | No stall |  -0.380|
| YDL060W   | No stall |  -0.382|
| YLR163C   | No stall |  -0.382|
| YDR212W   | Stall    |  -0.382|
| YNL102W   | Stall    |  -0.382|
| YOR367W   | No stall |  -0.382|
| YNL323W   | Stall    |  -0.383|
| YCR060W   | No stall |  -0.383|
| YGR204W   | No stall |  -0.383|
| YDR460W   | Stall    |  -0.384|
| YNL164C   | No stall |  -0.384|
| YMR117C   | No stall |  -0.384|
| YJL176C   | No stall |  -0.385|
| YKL215C   | No stall |  -0.385|
| YNL157W   | No stall |  -0.385|
| YHL009C   | No stall |  -0.385|
| YGR123C   | No stall |  -0.385|
| YHR046C   | No stall |  -0.385|
| YPL001W   | No stall |  -0.386|
| YPL215W   | No stall |  -0.386|
| YMR077C   | No stall |  -0.386|
| YML124C   | No stall |  -0.387|
| YIL140W   | No stall |  -0.388|
| YHL025W   | No stall |  -0.388|
| YOR311C   | No stall |  -0.389|
| YDR356W   | No stall |  -0.389|
| YBR265W   | No stall |  -0.389|
| YJR131W   | No stall |  -0.390|
| YLR363C   | No stall |  -0.390|
| YBR058C   | No stall |  -0.391|
| YIL007C   | No stall |  -0.391|
| YGR196C   | No stall |  -0.391|
| YLR056W   | No stall |  -0.391|
| YKR062W   | Stall    |  -0.392|
| YFR016C   | Stall    |  -0.392|
| YDR232W   | No stall |  -0.392|
| YLR066W   | No stall |  -0.392|
| YER073W   | No stall |  -0.393|
| YDR411C   | No stall |  -0.393|
| YMR006C   | No stall |  -0.393|
| YGL198W   | No stall |  -0.393|
| YJR064W   | No stall |  -0.395|
| YCR046C   | Stall    |  -0.395|
| YDR003W   | No stall |  -0.395|
| YPR190C   | No stall |  -0.396|
| YBR245C   | No stall |  -0.398|
| YNL183C   | No stall |  -0.398|
| YGR229C   | No stall |  -0.399|
| YDR101C   | No stall |  -0.399|
| YFL044C   | No stall |  -0.400|
| YMR060C   | No stall |  -0.400|
| YMR131C   | No stall |  -0.400|
| YDR006C   | No stall |  -0.400|
| YMR312W   | No stall |  -0.401|
| YKL195W   | No stall |  -0.401|
| YOR297C   | No stall |  -0.401|
| YLR197W   | Stall    |  -0.401|
| YOR112W   | No stall |  -0.401|
| YBR210W   | No stall |  -0.402|
| YDL173W   | Stall    |  -0.403|
| YBR286W   | No stall |  -0.403|
| YML015C   | No stall |  -0.404|
| YMR157C   | No stall |  -0.405|
| YOR281C   | No stall |  -0.405|
| YJL006C   | No stall |  -0.405|
| YMR093W   | Stall    |  -0.405|
| YJR068W   | No stall |  -0.406|
| YGR286C   | No stall |  -0.409|
| YHR135C   | No stall |  -0.409|
| YMR292W   | No stall |  -0.409|
| YGL129C   | No stall |  -0.409|
| YNL327W   | No stall |  -0.409|
| YKL114C   | Stall    |  -0.410|
| YIL103W   | No stall |  -0.411|
| YDR428C   | No stall |  -0.411|
| YDR244W   | No stall |  -0.412|
| YLR032W   | No stall |  -0.412|
| YGL246C   | No stall |  -0.412|
| YHR066W   | Stall    |  -0.412|
| YDL248W   | No stall |  -0.412|
| YGL120C   | Stall    |  -0.412|
| YDL205C   | No stall |  -0.412|
| YKL006C-A | No stall |  -0.413|
| YHR098C   | No stall |  -0.413|
| YER178W   | No stall |  -0.414|
| YGL012W   | No stall |  -0.414|
| YPL128C   | No stall |  -0.414|
| YGR048W   | Stall    |  -0.415|
| YDR449C   | Stall    |  -0.415|
| YDR470C   | No stall |  -0.415|
| YMR070W   | No stall |  -0.415|
| YBR162C   | No stall |  -0.416|
| YDL116W   | No stall |  -0.417|
| YPL096W   | No stall |  -0.417|
| YPR168W   | No stall |  -0.419|
| YLR274W   | No stall |  -0.419|
| YPL144W   | No stall |  -0.419|
| YJL186W   | No stall |  -0.420|
| YLL036C   | No stall |  -0.420|
| YNL100W   | No stall |  -0.420|
| YDR397C   | No stall |  -0.421|
| YGL137W   | No stall |  -0.421|
| YDR385W   | No stall |  -0.422|
| YBR236C   | No stall |  -0.422|
| YAL032C   | No stall |  -0.423|
| YBR283C   | No stall |  -0.423|
| YNL322C   | No stall |  -0.423|
| YIL044C   | No stall |  -0.424|
| YDL051W   | Stall    |  -0.424|
| YNL307C   | No stall |  -0.424|
| YOR158W   | No stall |  -0.424|
| YPR135W   | No stall |  -0.426|
| YBR145W   | No stall |  -0.426|
| YFL013C   | Stall    |  -0.426|
| YBR049C   | Stall    |  -0.426|
| YMR016C   | No stall |  -0.426|
| YMR125W   | No stall |  -0.426|
| YDR508C   | No stall |  -0.426|
| YBL009W   | No stall |  -0.426|
| YDL175C   | No stall |  -0.427|
| YBL018C   | No stall |  -0.427|
| YJR013W   | No stall |  -0.427|
| YAL059W   | Stall    |  -0.427|
| YDL212W   | No stall |  -0.427|
| YBL068W   | No stall |  -0.428|
| YDR117C   | Stall    |  -0.428|
| YDR354W   | No stall |  -0.429|
| YNR006W   | Stall    |  -0.430|
| YPR189W   | No stall |  -0.430|
| YJR057W   | No stall |  -0.431|
| YLR372W   | Stall    |  -0.432|
| YGL077C   | No stall |  -0.432|
| YER127W   | Stall    |  -0.432|
| YGL254W   | No stall |  -0.433|
| YPL030W   | No stall |  -0.433|
| YCL034W   | No stall |  -0.434|
| YML101C   | No stall |  -0.434|
| YLL034C   | Stall    |  -0.434|
| YPR024W   | No stall |  -0.434|
| YGR083C   | Stall    |  -0.434|
| YJL042W   | No stall |  -0.435|
| YGR245C   | Stall    |  -0.437|
| YNL233W   | Stall    |  -0.437|
| YPL118W   | No stall |  -0.437|
| YDL150W   | Stall    |  -0.437|
| YPL126W   | Stall    |  -0.437|
| YNL129W   | Stall    |  -0.437|
| YBR216C   | No stall |  -0.438|
| YJL033W   | Stall    |  -0.438|
| YGL203C   | Stall    |  -0.438|
| YHR181W   | No stall |  -0.439|
| YJL156C   | No stall |  -0.439|
| YDR196C   | No stall |  -0.439|
| YOR108W   | No stall |  -0.440|
| YKL020C   | Stall    |  -0.440|
| YER086W   | No stall |  -0.440|
| YNL154C   | No stall |  -0.441|
| YGR094W   | Stall    |  -0.442|
| YGR220C   | No stall |  -0.442|
| YMR139W   | No stall |  -0.444|
| YKL078W   | No stall |  -0.445|
| YLR351C   | No stall |  -0.446|
| YKR016W   | No stall |  -0.446|
| YJL050W   | Stall    |  -0.447|
| YDR152W   | No stall |  -0.448|
| YAL021C   | No stall |  -0.448|
| YOR119C   | Stall    |  -0.448|
| YGL061C   | No stall |  -0.448|
| YOL060C   | No stall |  -0.448|
| YMR115W   | No stall |  -0.449|
| YEL040W   | No stall |  -0.449|
| YOR165W   | No stall |  -0.449|
| YLR154C   | No stall |  -0.449|
| YNL123W   | No stall |  -0.449|
| YKL206C   | No stall |  -0.449|
| YBR211C   | No stall |  -0.450|
| YGR250C   | No stall |  -0.450|
| YDL143W   | No stall |  -0.451|
| YPL112C   | No stall |  -0.452|
| YLR201C   | No stall |  -0.452|
| YKL058W   | No stall |  -0.452|
| YOR133W   | No stall |  -0.453|
| YJL053W   | No stall |  -0.453|
| YNL191W   | No stall |  -0.454|
| YDR234W   | No stall |  -0.454|
| YER168C   | No stall |  -0.455|
| YKL125W   | No stall |  -0.455|
| YHR074W   | No stall |  -0.456|
| YDR080W   | Stall    |  -0.457|
| YKL166C   | No stall |  -0.458|
| YHR175W   | No stall |  -0.458|
| YDR448W   | No stall |  -0.459|
| YBL072C   | Stall    |  -0.460|
| YHR170W   | Stall    |  -0.460|
| YMR307W   | No stall |  -0.461|
| YDL230W   | No stall |  -0.462|
| YFL062W   | No stall |  -0.462|
| YJR043C   | No stall |  -0.462|
| YDR145W   | No stall |  -0.462|
| YOR039W   | No stall |  -0.463|
| YNL158W   | No stall |  -0.463|
| YGL190C   | No stall |  -0.463|
| YOL088C   | No stall |  -0.464|
| YJL130C   | No stall |  -0.464|
| YNL287W   | No stall |  -0.467|
| YBR095C   | No stall |  -0.467|
| YKR087C   | No stall |  -0.467|
| YHR114W   | No stall |  -0.467|
| YKL154W   | No stall |  -0.468|
| YGR246C   | Stall    |  -0.469|
| YDR430C   | No stall |  -0.469|
| YLR129W   | Stall    |  -0.469|
| YIL023C   | No stall |  -0.469|
| YNL066W   | No stall |  -0.469|
| YJR042W   | No stall |  -0.469|
| YDL132W   | No stall |  -0.470|
| YJR121W   | No stall |  -0.470|
| YMR186W   | Stall    |  -0.471|
| YLR245C   | No stall |  -0.471|
| YDL145C   | No stall |  -0.472|
| YOR250C   | No stall |  -0.473|
| YKL213C   | No stall |  -0.473|
| YKL140W   | No stall |  -0.473|
| YPL212C   | Stall    |  -0.473|
| YML115C   | No stall |  -0.476|
| YDL128W   | No stall |  -0.477|
| YGL164C   | No stall |  -0.477|
| YLR069C   | No stall |  -0.477|
| YOR195W   | Stall    |  -0.478|
| YEL065W   | No stall |  -0.479|
| YLR398C   | No stall |  -0.480|
| YOL152W   | No stall |  -0.480|
| YKL112W   | Stall    |  -0.480|
| YML097C   | No stall |  -0.480|
| YKL183W   | No stall |  -0.480|
| YBR289W   | Stall    |  -0.481|
| YGR177C   | No stall |  -0.481|
| YBR204C   | No stall |  -0.481|
| YPL259C   | Stall    |  -0.481|
| YGR101W   | No stall |  -0.481|
| YAL033W   | No stall |  -0.483|
| YDR380W   | No stall |  -0.483|
| YER047C   | Stall    |  -0.483|
| YDR297W   | No stall |  -0.484|
| YDL127W   | No stall |  -0.484|
| YBR093C   | No stall |  -0.484|
| YNR028W   | No stall |  -0.485|
| YCL045C   | No stall |  -0.487|
| YGL112C   | No stall |  -0.488|
| YHR005C   | No stall |  -0.488|
| YLR393W   | No stall |  -0.488|
| YNL306W   | No stall |  -0.488|
| YOR254C   | No stall |  -0.488|
| YHL030W   | No stall |  -0.489|
| YDR127W   | No stall |  -0.489|
| YJL074C   | No stall |  -0.490|
| YJL099W   | No stall |  -0.490|
| YOR155C   | No stall |  -0.490|
| YJL200C   | No stall |  -0.491|
| YDL179W   | No stall |  -0.492|
| YNL321W   | No stall |  -0.492|
| YGL053W   | No stall |  -0.492|
| YML127W   | No stall |  -0.492|
| YPL202C   | No stall |  -0.493|
| YHL007C   | Stall    |  -0.493|
| YLR373C   | No stall |  -0.493|
| YBR243C   | No stall |  -0.494|
| YFL023W   | Stall    |  -0.497|
| YML067C   | No stall |  -0.497|
| YFL017C   | No stall |  -0.497|
| YJR016C   | No stall |  -0.497|
| YJL062W   | No stall |  -0.497|
| YDL029W   | No stall |  -0.498|
| YHR034C   | No stall |  -0.498|
| YOR321W   | No stall |  -0.498|
| YDL102W   | Stall    |  -0.498|
| YNR023W   | No stall |  -0.498|
| YBR083W   | No stall |  -0.499|
| YPL023C   | No stall |  -0.499|
| YOR243C   | No stall |  -0.499|
| YOR262W   | No stall |  -0.499|
| YNL197C   | No stall |  -0.499|
| YDR324C   | Stall    |  -0.500|
| YBR238C   | No stall |  -0.500|
| YIL104C   | No stall |  -0.501|
| YNR067C   | No stall |  -0.503|
| YKL181W   | No stall |  -0.503|
| YDR046C   | No stall |  -0.503|
| YMR296C   | No stall |  -0.504|
| YKR090W   | No stall |  -0.505|
| YPR082C   | No stall |  -0.505|
| YBR092C   | No stall |  -0.507|
| YIL094C   | No stall |  -0.507|
| YLR440C   | No stall |  -0.507|
| YCL011C   | Stall    |  -0.508|
| YFL049W   | Stall    |  -0.508|
| YJR101W   | No stall |  -0.508|
| YJR051W   | No stall |  -0.509|
| YHR026W   | No stall |  -0.509|
| YBR169C   | No stall |  -0.510|
| YGR260W   | No stall |  -0.510|
| YGR165W   | No stall |  -0.511|
| YOR357C   | Stall    |  -0.512|
| YNL163C   | No stall |  -0.512|
| YBL035C   | No stall |  -0.512|
| YMR201C   | Stall    |  -0.513|
| YNR012W   | No stall |  -0.513|
| YDR231C   | No stall |  -0.514|
| YGL039W   | No stall |  -0.515|
| YDR027C   | No stall |  -0.517|
| YMR261C   | Stall    |  -0.517|
| YCL012C   | No stall |  -0.517|
| YFL048C   | No stall |  -0.517|
| YMR189W   | No stall |  -0.517|
| YOR083W   | No stall |  -0.518|
| YDR427W   | No stall |  -0.518|
| YDR539W   | No stall |  -0.518|
| YKL034W   | No stall |  -0.518|
| YDR329C   | No stall |  -0.519|
| YBL042C   | No stall |  -0.519|
| YDL020C   | No stall |  -0.519|
| YLR459W   | No stall |  -0.520|
| YHR133C   | No stall |  -0.520|
| YMR005W   | Stall    |  -0.520|
| YAL019W   | No stall |  -0.521|
| YMR250W   | No stall |  -0.521|
| YLR248W   | No stall |  -0.521|
| YCL027W   | Stall    |  -0.522|
| YIL111W   | No stall |  -0.522|
| YGL252C   | No stall |  -0.525|
| YLR226W   | No stall |  -0.525|
| YGL119W   | No stall |  -0.526|
| YPL100W   | No stall |  -0.526|
| YPL094C   | Stall    |  -0.526|
| YMR129W   | No stall |  -0.527|
| YGR231C   | No stall |  -0.528|
| YPL217C   | Stall    |  -0.528|
| YCL029C   | No stall |  -0.529|
| YMR149W   | No stall |  -0.530|
| YNL159C   | No stall |  -0.531|
| YFR033C   | No stall |  -0.531|
| YBL051C   | Stall    |  -0.531|
| YGR061C   | No stall |  -0.532|
| YKL080W   | No stall |  -0.532|
| YDR146C   | Stall    |  -0.533|
| YOR078W   | Stall    |  -0.534|
| YNL224C   | Stall    |  -0.535|
| YNL225C   | No stall |  -0.535|
| YDR398W   | Stall    |  -0.535|
| YOR057W   | Stall    |  -0.535|
| YPR120C   | No stall |  -0.536|
| YMR192W   | Stall    |  -0.536|
| YKL128C   | No stall |  -0.536|
| YPR023C   | No stall |  -0.537|
| YER053C   | No stall |  -0.538|
| YBR281C   | No stall |  -0.539|
| YER021W   | No stall |  -0.539|
| YOL031C   | No stall |  -0.539|
| YLR324W   | Stall    |  -0.539|
| YOR315W   | No stall |  -0.540|
| YER161C   | Stall    |  -0.541|
| YGL111W   | Stall    |  -0.541|
| YCL016C   | Stall    |  -0.541|
| YIL050W   | No stall |  -0.542|
| YDR005C   | No stall |  -0.542|
| YLR064W   | No stall |  -0.542|
| YAL051W   | No stall |  -0.542|
| YPR088C   | No stall |  -0.543|
| YHR101C   | No stall |  -0.545|
| YPR058W   | No stall |  -0.545|
| YJL128C   | No stall |  -0.545|
| YLR222C   | No stall |  -0.545|
| YDL052C   | No stall |  -0.546|
| YKL135C   | No stall |  -0.546|
| YML071C   | No stall |  -0.547|
| YJR144W   | Stall    |  -0.547|
| YML106W   | No stall |  -0.547|
| YGL095C   | No stall |  -0.547|
| YBR105C   | No stall |  -0.549|
| YLL006W   | No stall |  -0.549|
| YJL164C   | No stall |  -0.550|
| YOL002C   | Stall    |  -0.551|
| YCR030C   | Stall    |  -0.551|
| YEL042W   | No stall |  -0.551|
| YKL113C   | No stall |  -0.552|
| YGL038C   | No stall |  -0.552|
| YGR044C   | No stall |  -0.552|
| YGL065C   | No stall |  -0.553|
| YDR287W   | No stall |  -0.553|
| YER174C   | No stall |  -0.553|
| YOR125C   | No stall |  -0.553|
| YDR372C   | No stall |  -0.553|
| YBR212W   | No stall |  -0.553|
| YDR130C   | No stall |  -0.554|
| YLR357W   | Stall    |  -0.556|
| YLR079W   | No stall |  -0.556|
| YPR113W   | No stall |  -0.556|
| YNL148C   | No stall |  -0.556|
| YOL066C   | No stall |  -0.557|
| YGR138C   | Stall    |  -0.557|
| YHR032W   | No stall |  -0.558|
| YPR010C   | No stall |  -0.558|
| YBR193C   | No stall |  -0.558|
| YPL093W   | No stall |  -0.558|
| YAL023C   | No stall |  -0.559|
| YBR161W   | No stall |  -0.559|
| YER102W   | Stall    |  -0.560|
| YNL072W   | Stall    |  -0.562|
| YNR038W   | No stall |  -0.562|
| YGL125W   | No stall |  -0.562|
| YGL241W   | No stall |  -0.562|
| YPL087W   | No stall |  -0.563|
| YGL206C   | No stall |  -0.563|
| YDR169C   | No stall |  -0.564|
| YGR056W   | Stall    |  -0.564|
| YBR272C   | No stall |  -0.566|
| YPR162C   | No stall |  -0.567|
| YDR305C   | No stall |  -0.567|
| YJL061W   | No stall |  -0.568|
| YJL212C   | No stall |  -0.568|
| YMR235C   | No stall |  -0.569|
| YDR174W   | Stall    |  -0.571|
| YNL137C   | No stall |  -0.571|
| YHR198C   | No stall |  -0.572|
| YDL140C   | Stall    |  -0.572|
| YBL057C   | No stall |  -0.572|
| YMR145C   | No stall |  -0.573|
| YDR243C   | No stall |  -0.573|
| YOR267C   | Stall    |  -0.573|
| YGL223C   | No stall |  -0.573|
| YKL179C   | No stall |  -0.573|
| YLL021W   | No stall |  -0.574|
| YNL253W   | No stall |  -0.575|
| YHR056C   | No stall |  -0.575|
| YLR045C   | No stall |  -0.575|
| YHR024C   | No stall |  -0.576|
| YKL207W   | No stall |  -0.577|
| YLR378C   | No stall |  -0.577|
| YLL015W   | No stall |  -0.577|
| YDR331W   | No stall |  -0.577|
| YML125C   | No stall |  -0.578|
| YBR029C   | No stall |  -0.579|
| YJL192C   | No stall |  -0.579|
| YER183C   | No stall |  -0.580|
| YIL005W   | Stall    |  -0.580|
| YBR222C   | No stall |  -0.580|
| YFL036W   | No stall |  -0.581|
| YIL079C   | No stall |  -0.582|
| YBL037W   | No stall |  -0.582|
| YMR003W   | No stall |  -0.583|
| YER123W   | No stall |  -0.585|
| YNL090W   | No stall |  -0.585|
| YNR013C   | No stall |  -0.586|
| YBL102W   | No stall |  -0.586|
| YFR048W   | No stall |  -0.586|
| YLR342W   | Stall    |  -0.586|
| YDL195W   | No stall |  -0.586|
| YNR011C   | No stall |  -0.588|
| YDR292C   | Stall    |  -0.588|
| YJL157C   | No stall |  -0.588|
| YGR029W   | Stall    |  -0.589|
| YOR066W   | No stall |  -0.589|
| YBR068C   | No stall |  -0.590|
| YCL031C   | Stall    |  -0.590|
| YDL112W   | No stall |  -0.590|
| YOR023C   | Stall    |  -0.591|
| YMR012W   | No stall |  -0.592|
| YOL021C   | No stall |  -0.593|
| YDR312W   | Stall    |  -0.593|
| YMR297W   | No stall |  -0.594|
| YHR115C   | No stall |  -0.594|
| YJL141C   | No stall |  -0.596|
| YOR373W   | No stall |  -0.596|
| YCL057W   | Stall    |  -0.596|
| YFR042W   | No stall |  -0.598|
| YPL084W   | No stall |  -0.601|
| YJR082C   | No stall |  -0.602|
| YDR480W   | No stall |  -0.603|
| YPL049C   | No stall |  -0.603|
| YPL110C   | No stall |  -0.603|
| YDR072C   | Stall    |  -0.604|
| YER019W   | No stall |  -0.604|
| YDL042C   | No stall |  -0.605|
| YOR236W   | No stall |  -0.605|
| YJR075W   | No stall |  -0.605|
| YJR096W   | No stall |  -0.605|
| YGL021W   | No stall |  -0.606|
| YPL016W   | Stall    |  -0.606|
| YDL207W   | Stall    |  -0.606|
| YOL076W   | No stall |  -0.606|
| YKL094W   | No stall |  -0.606|
| YOR106W   | No stall |  -0.607|
| YMR275C   | Stall    |  -0.607|
| YPL101W   | No stall |  -0.607|
| YIL003W   | No stall |  -0.608|
| YPL018W   | No stall |  -0.609|
| YNL136W   | Stall    |  -0.610|
| YJR040W   | No stall |  -0.610|
| YJR080C   | No stall |  -0.611|
| YMR213W   | Stall    |  -0.611|
| YBR004C   | No stall |  -0.613|
| YLR214W   | No stall |  -0.613|
| YKR002W   | Stall    |  -0.614|
| YDR279W   | No stall |  -0.614|
| YPR017C   | No stall |  -0.614|
| YGR010W   | No stall |  -0.616|
| YJR132W   | No stall |  -0.616|
| YBR069C   | No stall |  -0.616|
| YOL122C   | No stall |  -0.618|
| YLR033W   | Stall    |  -0.618|
| YLR015W   | Stall    |  -0.619|
| YKL089W   | Stall    |  -0.619|
| YJR010W   | No stall |  -0.619|
| YKL208W   | No stall |  -0.620|
| YEL056W   | No stall |  -0.620|
| YBR079C   | Stall    |  -0.621|
| YHR169W   | Stall    |  -0.621|
| YPL053C   | No stall |  -0.622|
| YKL039W   | No stall |  -0.622|
| YDL063C   | No stall |  -0.622|
| YIL115C   | No stall |  -0.622|
| YOR132W   | No stall |  -0.623|
| YNL282W   | No stall |  -0.624|
| YGL169W   | No stall |  -0.624|
| YBR110W   | No stall |  -0.624|
| YNL315C   | No stall |  -0.625|
| YOR160W   | No stall |  -0.625|
| YDL134C   | No stall |  -0.626|
| YPR160W   | Stall    |  -0.626|
| YDL090C   | No stall |  -0.626|
| YIL110W   | No stall |  -0.629|
| YDR170C   | No stall |  -0.630|
| YDR507C   | Stall    |  -0.631|
| YPL183C   | No stall |  -0.631|
| YHR042W   | No stall |  -0.631|
| YJL004C   | No stall |  -0.631|
| YER055C   | No stall |  -0.631|
| YNL130C   | No stall |  -0.632|
| YPR042C   | Stall    |  -0.632|
| YAL041W   | No stall |  -0.633|
| YPR118W   | No stall |  -0.633|
| YER016W   | No stall |  -0.633|
| YKR037C   | No stall |  -0.634|
| YAL053W   | No stall |  -0.634|
| YOR174W   | No stall |  -0.636|
| YBR172C   | Stall    |  -0.636|
| YDL101C   | No stall |  -0.636|
| YER114C   | Stall    |  -0.636|
| YDR351W   | No stall |  -0.636|
| YLR254C   | No stall |  -0.637|
| YPR119W   | No stall |  -0.637|
| YNL039W   | Stall    |  -0.637|
| YMR237W   | No stall |  -0.637|
| YOR099W   | No stall |  -0.640|
| YDR295C   | Stall    |  -0.640|
| YKR085C   | Stall    |  -0.641|
| YDL008W   | No stall |  -0.642|
| YGR002C   | No stall |  -0.642|
| YLR249W   | Stall    |  -0.643|
| YGR200C   | No stall |  -0.644|
| YDL171C   | No stall |  -0.645|
| YBR302C   | No stall |  -0.645|
| YML132W   | No stall |  -0.645|
| YLR052W   | Stall    |  -0.645|
| YDR150W   | No stall |  -0.645|
| YER088C   | No stall |  -0.646|
| YML077W   | No stall |  -0.646|
| YAR008W   | No stall |  -0.647|
| YDR400W   | No stall |  -0.648|
| YBR123C   | No stall |  -0.648|
| YBR087W   | No stall |  -0.648|
| YJR112W   | No stall |  -0.650|
| YBL097W   | Stall    |  -0.650|
| YMR223W   | No stall |  -0.650|
| YBR278W   | No stall |  -0.650|
| YLR258W   | No stall |  -0.651|
| YKR088C   | No stall |  -0.651|
| YOR308C   | Stall    |  -0.651|
| YDR293C   | Stall    |  -0.651|
| YJL031C   | No stall |  -0.652|
| YDR060W   | Stall    |  -0.652|
| YEL055C   | Stall    |  -0.653|
| YLR397C   | No stall |  -0.653|
| YNL336W   | No stall |  -0.655|
| YAR014C   | No stall |  -0.656|
| YIL144W   | No stall |  -0.656|
| YBR017C   | No stall |  -0.656|
| YOR007C   | No stall |  -0.657|
| YDL098C   | No stall |  -0.657|
| YHR205W   | No stall |  -0.657|
| YJL065C   | No stall |  -0.658|
| YGR186W   | Stall    |  -0.658|
| YJL091C   | No stall |  -0.658|
| YIL155C   | No stall |  -0.658|
| YDR416W   | No stall |  -0.659|
| YNL201C   | No stall |  -0.660|
| YCR027C   | No stall |  -0.661|
| YOL159C   | No stall |  -0.661|
| YOR221C   | No stall |  -0.661|
| YML110C   | No stall |  -0.661|
| YHR108W   | No stall |  -0.662|
| YHR102W   | No stall |  -0.662|
| YJR033C   | No stall |  -0.662|
| YCR057C   | No stall |  -0.666|
| YER008C   | Stall    |  -0.666|
| YPR144C   | No stall |  -0.667|
| YOR138C   | Stall    |  -0.667|
| YDR245W   | No stall |  -0.669|
| YEL050C   | Stall    |  -0.669|
| YKR007W   | No stall |  -0.669|
| YOR110W   | No stall |  -0.670|
| YLL013C   | Stall    |  -0.670|
| YNL156C   | No stall |  -0.672|
| YMR001C   | No stall |  -0.672|
| YEL061C   | No stall |  -0.672|
| YMR267W   | No stall |  -0.673|
| YDR173C   | No stall |  -0.673|
| YBL022C   | No stall |  -0.673|
| YLR409C   | Stall    |  -0.673|
| YMR058W   | No stall |  -0.674|
| YBL089W   | No stall |  -0.675|
| YCR069W   | No stall |  -0.675|
| YFL025C   | No stall |  -0.676|
| YIL046W   | No stall |  -0.676|
| YPR148C   | No stall |  -0.677|
| YKL052C   | No stall |  -0.678|
| YJR006W   | No stall |  -0.681|
| YJR041C   | No stall |  -0.682|
| YOR075W   | No stall |  -0.682|
| YGL100W   | No stall |  -0.682|
| YPL009C   | Stall    |  -0.682|
| YML011C   | No stall |  -0.684|
| YKR072C   | No stall |  -0.684|
| YNL167C   | Stall    |  -0.684|
| YJR113C   | No stall |  -0.685|
| YOR150W   | No stall |  -0.685|
| YCR094W   | No stall |  -0.686|
| YBR002C   | No stall |  -0.687|
| YJL131C   | Stall    |  -0.688|
| YKL175W   | No stall |  -0.688|
| YLR059C   | No stall |  -0.689|
| YNR030W   | No stall |  -0.689|
| YLL004W   | No stall |  -0.690|
| YNL087W   | No stall |  -0.690|
| YGR143W   | No stall |  -0.691|
| YOR270C   | No stall |  -0.691|
| YFR021W   | No stall |  -0.691|
| YDR221W   | Stall    |  -0.692|
| YGR295C   | No stall |  -0.692|
| YMR228W   | No stall |  -0.692|
| YDR166C   | No stall |  -0.692|
| YBR086C   | No stall |  -0.692|
| YKR049C   | No stall |  -0.692|
| YGL110C   | Stall    |  -0.693|
| YAL026C   | No stall |  -0.694|
| YMR220W   | No stall |  -0.694|
| YDR501W   | Stall    |  -0.695|
| YGL107C   | No stall |  -0.695|
| YFL009W   | No stall |  -0.695|
| YBL006C   | Stall    |  -0.695|
| YJL002C   | No stall |  -0.696|
| YEL036C   | No stall |  -0.696|
| YOR115C   | No stall |  -0.696|
| YDL122W   | Stall    |  -0.698|
| YNL026W   | No stall |  -0.698|
| YLR083C   | No stall |  -0.698|
| YJL180C   | No stall |  -0.699|
| YHR036W   | No stall |  -0.700|
| YJL085W   | No stall |  -0.700|
| YLL011W   | Stall    |  -0.701|
| YKL218C   | No stall |  -0.701|
| YDR239C   | No stall |  -0.701|
| YMR061W   | No stall |  -0.701|
| YOR211C   | No stall |  -0.701|
| YGL207W   | No stall |  -0.702|
| YDR469W   | No stall |  -0.704|
| YML019W   | No stall |  -0.706|
| YNL124W   | Stall    |  -0.706|
| YNL005C   | No stall |  -0.706|
| YDR435C   | No stall |  -0.707|
| YDR268W   | No stall |  -0.707|
| YOR124C   | No stall |  -0.707|
| YDL188C   | No stall |  -0.708|
| YPR004C   | No stall |  -0.708|
| YNL286W   | No stall |  -0.709|
| YPR173C   | No stall |  -0.709|
| YJL082W   | No stall |  -0.709|
| YBR042C   | No stall |  -0.710|
| YFR025C   | No stall |  -0.710|
| YLR018C   | No stall |  -0.710|
| YJR059W   | No stall |  -0.711|
| YMR239C   | No stall |  -0.711|
| YHR144C   | No stall |  -0.711|
| YKR052C   | No stall |  -0.712|
| YDR505C   | No stall |  -0.713|
| YIL162W   | No stall |  -0.713|
| YDL031W   | Stall    |  -0.714|
| YKL059C   | No stall |  -0.715|
| YER125W   | No stall |  -0.715|
| YHR061C   | No stall |  -0.717|
| YHR007C   | No stall |  -0.717|
| YMR202W   | No stall |  -0.717|
| YHR062C   | No stall |  -0.718|
| YBL007C   | Stall    |  -0.719|
| YLR021W   | No stall |  -0.720|
| YGR013W   | No stall |  -0.720|
| YPR158W   | No stall |  -0.720|
| YBR132C   | No stall |  -0.721|
| YHR044C   | No stall |  -0.722|
| YMR037C   | Stall    |  -0.723|
| YPL210C   | Stall    |  -0.724|
| YDL113C   | No stall |  -0.725|
| YKR067W   | No stall |  -0.725|
| YLR386W   | No stall |  -0.725|
| YDL200C   | No stall |  -0.726|
| YDR517W   | Stall    |  -0.727|
| YBR262C   | No stall |  -0.727|
| YNR015W   | No stall |  -0.728|
| YER069W   | No stall |  -0.728|
| YKL073W   | Stall    |  -0.728|
| YNL280C   | No stall |  -0.728|
| YLR207W   | No stall |  -0.729|
| YCR048W   | No stall |  -0.729|
| YDR107C   | No stall |  -0.730|
| YNL219C   | No stall |  -0.730|
| YPL097W   | No stall |  -0.731|
| YBR200W   | No stall |  -0.731|
| YOR356W   | No stall |  -0.735|
| YLR291C   | No stall |  -0.735|
| YGR120C   | No stall |  -0.735|
| YNL048W   | No stall |  -0.735|
| YOL045W   | No stall |  -0.736|
| YDR165W   | No stall |  -0.736|
| YBR001C   | No stall |  -0.736|
| YBR202W   | No stall |  -0.737|
| YDR238C   | Stall    |  -0.738|
| YIL075C   | Stall    |  -0.738|
| YKR082W   | No stall |  -0.739|
| YGL027C   | No stall |  -0.740|
| YJL165C   | No stall |  -0.740|
| YMR247C   | No stall |  -0.740|
| YLR306W   | No stall |  -0.740|
| YIL048W   | No stall |  -0.741|
| YJR002W   | Stall    |  -0.741|
| YGL228W   | Stall    |  -0.743|
| YNL021W   | Stall    |  -0.744|
| YNL088W   | Stall    |  -0.745|
| YKL209C   | No stall |  -0.745|
| YGR227W   | Stall    |  -0.746|
| YIL043C   | No stall |  -0.746|
| YFL004W   | Stall    |  -0.746|
| YOR299W   | No stall |  -0.747|
| YPR128C   | No stall |  -0.747|
| YHR023W   | No stall |  -0.747|
| YKL165C   | No stall |  -0.748|
| YLR146C   | No stall |  -0.749|
| YDR062W   | No stall |  -0.749|
| YPL125W   | No stall |  -0.750|
| YIL114C   | No stall |  -0.751|
| YKL022C   | No stall |  -0.751|
| YER142C   | Stall    |  -0.751|
| YNR027W   | No stall |  -0.751|
| YML075C   | No stall |  -0.751|
| YHR030C   | No stall |  -0.752|
| YLR276C   | Stall    |  -0.753|
| YCR042C   | Stall    |  -0.753|
| YGL139W   | No stall |  -0.753|
| YBR115C   | No stall |  -0.753|
| YPL184C   | No stall |  -0.754|
| YGR102C   | No stall |  -0.754|
| YOR188W   | No stall |  -0.755|
| YOR037W   | No stall |  -0.755|
| YMR047C   | No stall |  -0.755|
| YGR252W   | No stall |  -0.755|
| YHR196W   | No stall |  -0.756|
| YJL019W   | No stall |  -0.756|
| YNL097C   | Stall    |  -0.756|
| YMR088C   | No stall |  -0.758|
| YLR277C   | No stall |  -0.760|
| YER040W   | Stall    |  -0.760|
| YDR147W   | No stall |  -0.762|
| YOR076C   | No stall |  -0.763|
| YBL098W   | No stall |  -0.764|
| YML121W   | No stall |  -0.765|
| YOR087W   | Stall    |  -0.765|
| YPL254W   | No stall |  -0.766|
| YBR201W   | No stall |  -0.766|
| YBR146W   | Stall    |  -0.767|
| YIR010W   | Stall    |  -0.768|
| YLR078C   | No stall |  -0.768|
| YER014W   | No stall |  -0.768|
| YLR182W   | No stall |  -0.769|
| YPL138C   | Stall    |  -0.769|
| YLL029W   | Stall    |  -0.770|
| YER116C   | Stall    |  -0.771|
| YLR219W   | Stall    |  -0.771|
| YGR012W   | No stall |  -0.773|
| YDL088C   | No stall |  -0.773|
| YDR294C   | No stall |  -0.774|
| YPL256C   | No stall |  -0.774|
| YLR215C   | No stall |  -0.774|
| YOR301W   | No stall |  -0.777|
| YHR050W   | No stall |  -0.777|
| YER001W   | No stall |  -0.778|
| YBR199W   | No stall |  -0.779|
| YOL090W   | Stall    |  -0.779|
| YMR319C   | No stall |  -0.779|
| YDL203C   | No stall |  -0.780|
| YPR171W   | Stall    |  -0.780|
| YGL093W   | No stall |  -0.780|
| YJR161C   | No stall |  -0.781|
| YPR056W   | Stall    |  -0.781|
| YLR212C   | No stall |  -0.782|
| YKL182W   | No stall |  -0.782|
| YGR233C   | No stall |  -0.782|
| YPR185W   | No stall |  -0.783|
| YER151C   | Stall    |  -0.783|
| YJR049C   | No stall |  -0.784|
| YHR197W   | No stall |  -0.785|
| YFR007W   | No stall |  -0.785|
| YJR058C   | No stall |  -0.785|
| YPL038W   | No stall |  -0.786|
| YGL151W   | No stall |  -0.786|
| YGR103W   | Stall    |  -0.786|
| YPR075C   | Stall    |  -0.787|
| YPR127W   | No stall |  -0.787|
| YFL026W   | No stall |  -0.788|
| YJL168C   | Stall    |  -0.788|
| YBR151W   | No stall |  -0.788|
| YNL240C   | No stall |  -0.790|
| YDR110W   | Stall    |  -0.790|
| YOR260W   | No stall |  -0.790|
| YOR219C   | No stall |  -0.790|
| YOL042W   | No stall |  -0.791|
| YGR100W   | No stall |  -0.791|
| YLL061W   | No stall |  -0.793|
| YHR003C   | Stall    |  -0.793|
| YDR184C   | No stall |  -0.794|
| YPL241C   | No stall |  -0.794|
| YMR029C   | No stall |  -0.794|
| YNL085W   | No stall |  -0.795|
| YFL041W   | No stall |  -0.795|
| YGL227W   | No stall |  -0.795|
| YPL019C   | Stall    |  -0.798|
| YMR140W   | No stall |  -0.798|
| YGL150C   | Stall    |  -0.798|
| YBR103W   | No stall |  -0.798|
| YOR207C   | No stall |  -0.799|
| YML111W   | Stall    |  -0.800|
| YDL111C   | No stall |  -0.801|
| YGL196W   | No stall |  -0.802|
| YBL074C   | No stall |  -0.802|
| YDL093W   | No stall |  -0.803|
| YOR372C   | Stall    |  -0.805|
| YLR096W   | No stall |  -0.805|
| YAL042W   | No stall |  -0.806|
| YER149C   | No stall |  -0.806|
| YDR335W   | No stall |  -0.806|
| YLR131C   | No stall |  -0.807|
| YJL081C   | No stall |  -0.807|
| YDR464W   | Stall    |  -0.807|
| YDR013W   | No stall |  -0.807|
| YJL201W   | No stall |  -0.808|
| YOR098C   | No stall |  -0.809|
| YKR038C   | No stall |  -0.810|
| YHR161C   | No stall |  -0.811|
| YER105C   | No stall |  -0.812|
| YMR086W   | Stall    |  -0.814|
| YOR217W   | Stall    |  -0.814|
| YOR130C   | No stall |  -0.814|
| YBR070C   | No stall |  -0.815|
| YMR165C   | No stall |  -0.817|
| YCR065W   | No stall |  -0.817|
| YPR029C   | No stall |  -0.817|
| YER078C   | No stall |  -0.820|
| YPR175W   | No stall |  -0.820|
| YNL049C   | No stall |  -0.821|
| YDR531W   | No stall |  -0.821|
| YGR247W   | No stall |  -0.822|
| YNL051W   | No stall |  -0.822|
| YHR186C   | No stall |  -0.823|
| YPR048W   | No stall |  -0.823|
| YDR311W   | No stall |  -0.824|
| YPL226W   | Stall    |  -0.825|
| YNL041C   | No stall |  -0.825|
| YJR090C   | No stall |  -0.825|
| YDR303C   | No stall |  -0.826|
| YLR113W   | No stall |  -0.827|
| YKL079W   | No stall |  -0.827|
| YJR005W   | No stall |  -0.828|
| YJL198W   | No stall |  -0.828|
| YPL273W   | No stall |  -0.828|
| YJR100C   | No stall |  -0.828|
| YLR090W   | No stall |  -0.829|
| YIL009C-A | No stall |  -0.830|
| YOR153W   | Stall    |  -0.831|
| YLL012W   | No stall |  -0.831|
| YOR241W   | No stall |  -0.832|
| YJL162C   | Stall    |  -0.832|
| YGL017W   | No stall |  -0.833|
| YKL146W   | No stall |  -0.833|
| YMR289W   | No stall |  -0.834|
| YBR157C   | No stall |  -0.835|
| YLR166C   | No stall |  -0.837|
| YPL057C   | No stall |  -0.838|
| YBR023C   | No stall |  -0.838|
| YJR001W   | No stall |  -0.839|
| YPL116W   | No stall |  -0.839|
| YMR113W   | No stall |  -0.840|
| YNL099C   | No stall |  -0.840|
| YGL163C   | No stall |  -0.840|
| YIL119C   | Stall    |  -0.841|
| YOR048C   | Stall    |  -0.842|
| YKL149C   | No stall |  -0.843|
| YDR181C   | Stall    |  -0.843|
| YGR091W   | Stall    |  -0.843|
| YGL195W   | No stall |  -0.845|
| YMR243C   | No stall |  -0.846|
| YKL174C   | No stall |  -0.847|
| YDL224C   | No stall |  -0.849|
| YGL083W   | No stall |  -0.850|
| YDR034C   | Stall    |  -0.850|
| YCR082W   | No stall |  -0.853|
| YPL224C   | No stall |  -0.853|
| YNL258C   | No stall |  -0.853|
| YBR084W   | No stall |  -0.854|
| YNL169C   | Stall    |  -0.854|
| YDL005C   | No stall |  -0.856|
| YLR330W   | Stall    |  -0.856|
| YOL051W   | No stall |  -0.856|
| YHR163W   | No stall |  -0.856|
| YGL201C   | No stall |  -0.856|
| YDR105C   | No stall |  -0.857|
| YOR233W   | Stall    |  -0.858|
| YDL159W   | Stall    |  -0.859|
| YOL145C   | Stall    |  -0.860|
| YJL098W   | No stall |  -0.860|
| YNL182C   | No stall |  -0.860|
| YHL020C   | No stall |  -0.861|
| YGL086W   | No stall |  -0.862|
| YDL141W   | No stall |  -0.862|
| YIL030C   | No stall |  -0.862|
| YPL050C   | No stall |  -0.863|
| YKR100C   | No stall |  -0.863|
| YDR236C   | No stall |  -0.864|
| YBR150C   | Stall    |  -0.866|
| YGR178C   | No stall |  -0.866|
| YBR073W   | No stall |  -0.867|
| YNL298W   | No stall |  -0.869|
| YDR204W   | Stall    |  -0.869|
| YAL043C   | No stall |  -0.870|
| YOR360C   | No stall |  -0.870|
| YBL101C   | No stall |  -0.871|
| YOL089C   | Stall    |  -0.871|
| YBR015C   | Stall    |  -0.871|
| YGL224C   | No stall |  -0.872|
| YML109W   | Stall    |  -0.874|
| YDR007W   | No stall |  -0.874|
| YNL006W   | No stall |  -0.874|
| YKL116C   | No stall |  -0.875|
| YGR271W   | No stall |  -0.876|
| YPL032C   | Stall    |  -0.877|
| YOR304W   | No stall |  -0.877|
| YMR032W   | No stall |  -0.877|
| YNL221C   | Stall    |  -0.878|
| YDR410C   | No stall |  -0.879|
| YPL246C   | No stall |  -0.879|
| YLR088W   | No stall |  -0.880|
| YGL173C   | Stall    |  -0.880|
| YDR103W   | No stall |  -0.880|
| YLL035W   | No stall |  -0.881|
| YGR104C   | No stall |  -0.882|
| YCR044C   | No stall |  -0.883|
| YDR143C   | No stall |  -0.884|
| YDR208W   | Stall    |  -0.884|
| YHR172W   | No stall |  -0.884|
| YOR353C   | No stall |  -0.884|
| YIL088C   | No stall |  -0.884|
| YBL082C   | No stall |  -0.884|
| YGL160W   | No stall |  -0.885|
| YBR130C   | No stall |  -0.886|
| YBL055C   | No stall |  -0.887|
| YMR308C   | No stall |  -0.888|
| YER013W   | Stall    |  -0.889|
| YJL139C   | No stall |  -0.889|
| YDR479C   | No stall |  -0.889|
| YCR017C   | No stall |  -0.891|
| YBR142W   | Stall    |  -0.894|
| YNR074C   | No stall |  -0.894|
| YLR206W   | Stall    |  -0.894|
| YDR135C   | No stall |  -0.895|
| YJR024C   | No stall |  -0.895|
| YMR073C   | No stall |  -0.895|
| YDR434W   | No stall |  -0.896|
| YCR028C   | No stall |  -0.896|
| YIL068C   | No stall |  -0.896|
| YBR208C   | Stall    |  -0.897|
| YFL047W   | No stall |  -0.898|
| YOR175C   | No stall |  -0.898|
| YKL057C   | No stall |  -0.899|
| YLR006C   | No stall |  -0.899|
| YOL003C   | No stall |  -0.899|
| YDR257C   | No stall |  -0.901|
| YNL256W   | No stall |  -0.903|
| YDR376W   | No stall |  -0.905|
| YBL085W   | No stall |  -0.905|
| YMR164C   | Stall    |  -0.906|
| YML031W   | No stall |  -0.907|
| YLR405W   | No stall |  -0.909|
| YOR081C   | No stall |  -0.910|
| YKL106W   | No stall |  -0.910|
| YOR229W   | No stall |  -0.912|
| YKL010C   | Stall    |  -0.912|
| YDR141C   | Stall    |  -0.914|
| YBR125C   | No stall |  -0.915|
| YNL059C   | Stall    |  -0.915|
| YHR206W   | No stall |  -0.916|
| YPL103C   | No stall |  -0.916|
| YJR143C   | No stall |  -0.917|
| YBR041W   | No stall |  -0.917|
| YOR090C   | No stall |  -0.917|
| YFR002W   | No stall |  -0.917|
| YDR017C   | No stall |  -0.918|
| YER171W   | No stall |  -0.918|
| YMR273C   | Stall    |  -0.919|
| YKR076W   | No stall |  -0.920|
| YKL129C   | Stall    |  -0.920|
| YPR106W   | No stall |  -0.920|
| YOR035C   | No stall |  -0.921|
| YGR007W   | No stall |  -0.921|
| YPL233W   | No stall |  -0.922|
| YDR530C   | No stall |  -0.924|
| YLR190W   | No stall |  -0.924|
| YBL079W   | No stall |  -0.924|
| YJR066W   | No stall |  -0.925|
| YFR038W   | Stall    |  -0.926|
| YLR071C   | No stall |  -0.926|
| YDR498C   | No stall |  -0.926|
| YDR367W   | No stall |  -0.927|
| YLR116W   | No stall |  -0.928|
| YHR072W   | No stall |  -0.928|
| YOR151C   | No stall |  -0.928|
| YFR005C   | No stall |  -0.928|
| YLR399C   | Stall    |  -0.930|
| YGR014W   | No stall |  -0.930|
| YGL067W   | No stall |  -0.930|
| YML056C   | No stall |  -0.930|
| YNR051C   | No stall |  -0.931|
| YDR211W   | Stall    |  -0.932|
| YBR260C   | No stall |  -0.933|
| YJL154C   | No stall |  -0.933|
| YEL013W   | No stall |  -0.933|
| YEL031W   | No stall |  -0.934|
| YOL116W   | No stall |  -0.936|
| YPL242C   | No stall |  -0.936|
| YMR229C   | Stall    |  -0.936|
| YPR179C   | Stall    |  -0.937|
| YPL180W   | No stall |  -0.937|
| YNL118C   | No stall |  -0.938|
| YJL035C   | No stall |  -0.938|
| YLR319C   | No stall |  -0.941|
| YMR266W   | Stall    |  -0.942|
| YHR039C   | No stall |  -0.942|
| YIL131C   | No stall |  -0.942|
| YBR175W   | No stall |  -0.942|
| YNL139C   | No stall |  -0.942|
| YBR205W   | No stall |  -0.943|
| YHL048W   | No stall |  -0.945|
| YOR067C   | No stall |  -0.946|
| YOR025W   | No stall |  -0.948|
| YGR274C   | Stall    |  -0.948|
| YLR410W   | Stall    |  -0.948|
| YER010C   | No stall |  -0.949|
| YFL021W   | Stall    |  -0.949|
| YCL061C   | Stall    |  -0.950|
| YCR011C   | No stall |  -0.950|
| YGR147C   | Stall    |  -0.952|
| YAL022C   | No stall |  -0.953|
| YNL299W   | Stall    |  -0.954|
| YBR059C   | No stall |  -0.954|
| YML005W   | No stall |  -0.954|
| YJL183W   | No stall |  -0.955|
| YJL039C   | No stall |  -0.956|
| YLR077W   | No stall |  -0.956|
| YAR042W   | No stall |  -0.956|
| YMR264W   | No stall |  -0.957|
| YGR041W   | No stall |  -0.958|
| YGR032W   | Stall    |  -0.958|
| YPL092W   | No stall |  -0.958|
| YHR083W   | No stall |  -0.959|
| YDR144C   | No stall |  -0.959|
| YFR028C   | Stall    |  -0.960|
| YCL052C   | Stall    |  -0.961|
| YNL267W   | No stall |  -0.963|
| YGR198W   | No stall |  -0.963|
| YPR186C   | Stall    |  -0.964|
| YMR219W   | Stall    |  -0.965|
| YMR277W   | No stall |  -0.966|
| YBR288C   | No stall |  -0.966|
| YDR441C   | No stall |  -0.969|
| YHR028C   | No stall |  -0.970|
| YOR278W   | No stall |  -0.970|
| YIR006C   | Stall    |  -0.971|
| YKL014C   | No stall |  -0.971|
| YLR085C   | No stall |  -0.972|
| YLR039C   | No stall |  -0.973|
| YDR140W   | No stall |  -0.973|
| YDL174C   | Stall    |  -0.974|
| YGR122W   | No stall |  -0.974|
| YIL095W   | Stall    |  -0.975|
| YEL048C   | No stall |  -0.975|
| YNL293W   | Stall    |  -0.976|
| YER164W   | Stall    |  -0.977|
| YMR224C   | Stall    |  -0.977|
| YGR281W   | Stall    |  -0.977|
| YDR475C   | Stall    |  -0.977|
| YMR304W   | No stall |  -0.978|
| YPL012W   | Stall    |  -0.979|
| YPR104C   | Stall    |  -0.979|
| YGR068C   | No stall |  -0.982|
| YIR018W   | No stall |  -0.982|
| YGL092W   | No stall |  -0.982|
| YDR261C   | No stall |  -0.983|
| YLL040C   | Stall    |  -0.983|
| YBL056W   | No stall |  -0.983|
| YDR207C   | Stall    |  -0.984|
| YPL105C   | Stall    |  -0.984|
| YJL034W   | No stall |  -0.984|
| YLR240W   | No stall |  -0.985|
| YPL152W   | No stall |  -0.985|
| YLR348C   | No stall |  -0.986|
| YLR095C   | Stall    |  -0.986|
| YLR107W   | No stall |  -0.987|
| YPR055W   | No stall |  -0.988|
| YPL227C   | No stall |  -0.989|
| YBR256C   | No stall |  -0.990|
| YLR016C   | No stall |  -0.991|
| YBR179C   | Stall    |  -0.992|
| YGR113W   | No stall |  -0.993|
| YMR199W   | No stall |  -0.994|
| YFR031C   | No stall |  -0.994|
| YKL019W   | No stall |  -0.996|
| YAL017W   | Stall    |  -0.996|
| YKL068W   | No stall |  -0.997|
| YLR117C   | Stall    |  -0.998|
| YGR046W   | No stall |  -1.000|
| YFR030W   | No stall |  -1.002|
| YML049C   | No stall |  -1.003|
| YGL219C   | Stall    |  -1.004|
| YDL058W   | No stall |  -1.006|
| YLR223C   | Stall    |  -1.006|
| YDR490C   | No stall |  -1.006|
| YHR188C   | No stall |  -1.006|
| YFR019W   | Stall    |  -1.007|
| YPL020C   | Stall    |  -1.009|
| YLL043W   | Stall    |  -1.009|
| YJR140C   | No stall |  -1.010|
| YGL005C   | No stall |  -1.011|
| YNL152W   | No stall |  -1.012|
| YOR377W   | No stall |  -1.012|
| YFL002C   | Stall    |  -1.012|
| YPL007C   | No stall |  -1.012|
| YNL077W   | Stall    |  -1.013|
| YJL110C   | No stall |  -1.013|
| YER075C   | No stall |  -1.014|
| YOR071C   | No stall |  -1.014|
| YCL063W   | No stall |  -1.014|
| YLR019W   | No stall |  -1.015|
| YDR235W   | Stall    |  -1.015|
| YAR018C   | No stall |  -1.017|
| YML059C   | No stall |  -1.017|
| YMR076C   | Stall    |  -1.017|
| YMR284W   | Stall    |  -1.019|
| YGL142C   | No stall |  -1.020|
| YNL227C   | Stall    |  -1.022|
| YBL076C   | Stall    |  -1.022|
| YGR255C   | No stall |  -1.022|
| YGL115W   | No stall |  -1.023|
| YJR062C   | No stall |  -1.023|
| YDR365C   | Stall    |  -1.024|
| YMR097C   | No stall |  -1.026|
| YFR040W   | Stall    |  -1.028|
| YDL237W   | No stall |  -1.029|
| YHR154W   | No stall |  -1.029|
| YDR409W   | No stall |  -1.029|
| YER089C   | No stall |  -1.029|
| YPR134W   | No stall |  -1.030|
| YNL180C   | Stall    |  -1.032|
| YCR019W   | No stall |  -1.032|
| YHL039W   | No stall |  -1.033|
| YNL326C   | No stall |  -1.033|
| YOR120W   | No stall |  -1.034|
| YOR256C   | No stall |  -1.035|
| YGR099W   | No stall |  -1.036|
| YLR403W   | No stall |  -1.037|
| YLR452C   | No stall |  -1.038|
| YGR003W   | No stall |  -1.039|
| YJL101C   | No stall |  -1.039|
| YKR003W   | No stall |  -1.040|
| YNL250W   | No stall |  -1.041|
| YJL207C   | Stall    |  -1.041|
| YPR169W   | Stall    |  -1.042|
| YOR147W   | No stall |  -1.043|
| YGL073W   | No stall |  -1.044|
| YLL048C   | No stall |  -1.044|
| YGR146C   | No stall |  -1.045|
| YLR450W   | No stall |  -1.047|
| YOR201C   | No stall |  -1.047|
| YMR231W   | No stall |  -1.047|
| YML128C   | No stall |  -1.047|
| YOR047C   | No stall |  -1.049|
| YGL248W   | No stall |  -1.052|
| YDR538W   | No stall |  -1.053|
| YNL333W   | No stall |  -1.055|
| YNL068C   | No stall |  -1.055|
| YDL025C   | Stall    |  -1.056|
| YNL078W   | Stall    |  -1.056|
| YDL189W   | Stall    |  -1.057|
| YGL062W   | Stall    |  -1.057|
| YIL153W   | No stall |  -1.058|
| YFR029W   | Stall    |  -1.058|
| YLR389C   | No stall |  -1.058|
| YMR080C   | No stall |  -1.058|
| YHR031C   | No stall |  -1.059|
| YNL119W   | No stall |  -1.060|
| YOL158C   | No stall |  -1.061|
| YHL040C   | No stall |  -1.062|
| YER166W   | Stall    |  -1.063|
| YGR040W   | No stall |  -1.063|
| YDR264C   | No stall |  -1.063|
| YER018C   | No stall |  -1.064|
| YDR439W   | Stall    |  -1.064|
| YIL107C   | Stall    |  -1.064|
| YNL103W   | Stall    |  -1.065|
| YLR382C   | No stall |  -1.068|
| YGR065C   | No stall |  -1.068|
| YPR140W   | No stall |  -1.069|
| YGR033C   | No stall |  -1.070|
| YOL020W   | No stall |  -1.070|
| YKL134C   | No stall |  -1.070|
| YIR023W   | No stall |  -1.073|
| YBR003W   | No stall |  -1.073|
| YER061C   | No stall |  -1.074|
| YNL199C   | No stall |  -1.074|
| YKL108W   | Stall    |  -1.074|
| YJL194W   | Stall    |  -1.075|
| YNR033W   | No stall |  -1.075|
| YIL121W   | No stall |  -1.076|
| YER109C   | No stall |  -1.076|
| YGL143C   | Stall    |  -1.076|
| YOR320C   | No stall |  -1.076|
| YOR086C   | Stall    |  -1.077|
| YBL061C   | Stall    |  -1.078|
| YML061C   | Stall    |  -1.082|
| YJL005W   | Stall    |  -1.082|
| YDR457W   | Stall    |  -1.084|
| YNL206C   | No stall |  -1.085|
| YKL173W   | Stall    |  -1.086|
| YBR253W   | No stall |  -1.087|
| YCL038C   | No stall |  -1.088|
| YIL128W   | No stall |  -1.088|
| YMR054W   | No stall |  -1.089|
| YMR302C   | No stall |  -1.090|
| YOR058C   | No stall |  -1.091|
| YLR238W   | Stall    |  -1.091|
| YOR322C   | No stall |  -1.092|
| YDL240W   | No stall |  -1.092|
| YGR108W   | No stall |  -1.093|
| YML103C   | No stall |  -1.093|
| YEL022W   | No stall |  -1.098|
| YJR091C   | No stall |  -1.098|
| YPL176C   | Stall    |  -1.100|
| YDL155W   | No stall |  -1.101|
| YDR456W   | No stall |  -1.102|
| YDL056W   | Stall    |  -1.102|
| YGL016W   | No stall |  -1.102|
| YMR089C   | Stall    |  -1.102|
| YJR031C   | No stall |  -1.103|
| YOR116C   | No stall |  -1.103|
| YFL059W   | No stall |  -1.107|
| YNL008C   | Stall    |  -1.108|
| YFL008W   | Stall    |  -1.110|
| YLR371W   | No stall |  -1.111|
| YNL297C   | No stall |  -1.112|
| YHR017W   | No stall |  -1.112|
| YDR497C   | No stall |  -1.114|
| YOR166C   | No stall |  -1.115|
| YOR001W   | Stall    |  -1.115|
| YDR093W   | Stall    |  -1.116|
| YPL270W   | No stall |  -1.118|
| YEL053C   | No stall |  -1.119|
| YGR028W   | No stall |  -1.121|
| YPL069C   | No stall |  -1.121|
| YHR199C   | No stall |  -1.121|
| YMR150C   | No stall |  -1.122|
| YBL017C   | Stall    |  -1.122|
| YPL002C   | No stall |  -1.123|
| YDR081C   | Stall    |  -1.125|
| YNR032W   | No stall |  -1.125|
| YDR202C   | No stall |  -1.125|
| YDR164C   | Stall    |  -1.125|
| YKL124W   | Stall    |  -1.126|
| YHR041C   | No stall |  -1.127|
| YOR371C   | No stall |  -1.127|
| YDR228C   | No stall |  -1.128|
| YDR414C   | No stall |  -1.128|
| YKR089C   | No stall |  -1.128|
| YOR290C   | Stall    |  -1.131|
| YKL012W   | No stall |  -1.131|
| YDL074C   | No stall |  -1.132|
| YOL049W   | No stall |  -1.133|
| YCR081W   | No stall |  -1.133|
| YHR086W   | No stall |  -1.134|
| YGR241C   | No stall |  -1.136|
| YBR129C   | Stall    |  -1.138|
| YHR110W   | No stall |  -1.142|
| YJR126C   | No stall |  -1.142|
| YLL031C   | No stall |  -1.143|
| YJL133W   | No stall |  -1.146|
| YIL112W   | Stall    |  -1.148|
| YPL045W   | No stall |  -1.148|
| YIR033W   | Stall    |  -1.151|
| YDR363W   | No stall |  -1.151|
| YGL014W   | Stall    |  -1.154|
| YKL191W   | No stall |  -1.154|
| YKR095W   | No stall |  -1.155|
| YDR122W   | Stall    |  -1.155|
| YBR081C   | Stall    |  -1.156|
| YDR347W   | No stall |  -1.156|
| YDL117W   | Stall    |  -1.156|
| YDR452W   | Stall    |  -1.157|
| YOR319W   | No stall |  -1.158|
| YDR375C   | No stall |  -1.158|
| YDR422C   | No stall |  -1.159|
| YBL020W   | No stall |  -1.160|
| YDR159W   | No stall |  -1.162|
| YEL029C   | No stall |  -1.163|
| YDL028C   | No stall |  -1.166|
| YKL198C   | Stall    |  -1.166|
| YBR097W   | No stall |  -1.170|
| YKL004W   | No stall |  -1.172|
| YOR144C   | No stall |  -1.172|
| YPL058C   | No stall |  -1.173|
| YDR194C   | No stall |  -1.173|
| YIL090W   | No stall |  -1.174|
| YJL073W   | No stall |  -1.174|
| YDR483W   | No stall |  -1.175|
| YAL048C   | Stall    |  -1.176|
| YCR068W   | No stall |  -1.176|
| YCL014W   | No stall |  -1.177|
| YEL032W   | Stall    |  -1.179|
| YKR019C   | Stall    |  -1.183|
| YPL083C   | No stall |  -1.183|
| YDR096W   | Stall    |  -1.183|
| YIR021W   | No stall |  -1.184|
| YNL003C   | No stall |  -1.185|
| YER172C   | No stall |  -1.185|
| YNL076W   | No stall |  -1.185|
| YKL205W   | No stall |  -1.186|
| YJL029C   | No stall |  -1.187|
| YLR237W   | No stall |  -1.187|
| YBR229C   | No stall |  -1.188|
| YHR195W   | No stall |  -1.188|
| YGL022W   | No stall |  -1.188|
| YKL028W   | No stall |  -1.189|
| YDL167C   | No stall |  -1.192|
| YAL001C   | Stall    |  -1.192|
| YNL275W   | No stall |  -1.192|
| YGR097W   | No stall |  -1.193|
| YJL112W   | Stall    |  -1.193|
| YGL023C   | No stall |  -1.194|
| YDR219C   | No stall |  -1.194|
| YLR384C   | Stall    |  -1.196|
| YGL238W   | No stall |  -1.197|
| YMR212C   | No stall |  -1.197|
| YNR031C   | No stall |  -1.198|
| YGR270W   | Stall    |  -1.199|
| YCL036W   | No stall |  -1.199|
| YCR059C   | No stall |  -1.201|
| YKR063C   | Stall    |  -1.208|
| YJL051W   | No stall |  -1.209|
| YGL211W   | No stall |  -1.211|
| YIL158W   | No stall |  -1.211|
| YML029W   | No stall |  -1.212|
| YMR013C   | No stall |  -1.212|
| YGL060W   | No stall |  -1.212|
| YOR092W   | No stall |  -1.213|
| YIL129C   | No stall |  -1.215|
| YJL094C   | No stall |  -1.215|
| YBL067C   | Stall    |  -1.219|
| YER157W   | No stall |  -1.220|
| YLR106C   | Stall    |  -1.222|
| YHR084W   | No stall |  -1.223|
| YBR218C   | Stall    |  -1.223|
| YER154W   | No stall |  -1.226|
| YML032C   | No stall |  -1.227|
| YLR105C   | No stall |  -1.227|
| YER017C   | No stall |  -1.228|
| YGR140W   | Stall    |  -1.229|
| YOR341W   | Stall    |  -1.229|
| YOR109W   | Stall    |  -1.233|
| YNL014W   | Stall    |  -1.234|
| YMR234W   | No stall |  -1.234|
| YDL070W   | Stall    |  -1.236|
| YPL151C   | No stall |  -1.237|
| YPL195W   | Stall    |  -1.237|
| YER113C   | No stall |  -1.238|
| YGR166W   | No stall |  -1.238|
| YIL014W   | No stall |  -1.240|
| YLR361C   | No stall |  -1.240|
| YPL042C   | No stall |  -1.242|
| YML116W   | No stall |  -1.243|
| YDL089W   | Stall    |  -1.244|
| YDL019C   | Stall    |  -1.249|
| YOL148C   | Stall    |  -1.251|
| YFL024C   | Stall    |  -1.251|
| YML114C   | No stall |  -1.253|
| YPL172C   | No stall |  -1.253|
| YNR047W   | Stall    |  -1.254|
| YMR216C   | Stall    |  -1.256|
| YJR110W   | No stall |  -1.256|
| YBR065C   | Stall    |  -1.257|
| YML072C   | Stall    |  -1.257|
| YCR075C   | No stall |  -1.257|
| YNL083W   | No stall |  -1.261|
| YMR035W   | No stall |  -1.261|
| YBL084C   | No stall |  -1.263|
| YNL262W   | Stall    |  -1.263|
| YGL002W   | No stall |  -1.263|
| YGL071W   | Stall    |  -1.263|
| YOR307C   | No stall |  -1.265|
| YDR301W   | Stall    |  -1.266|
| YGL216W   | No stall |  -1.268|
| YDL077C   | No stall |  -1.270|
| YNL127W   | No stall |  -1.271|
| YPR070W   | No stall |  -1.271|
| YKL155C   | Stall    |  -1.272|
| YLR228C   | No stall |  -1.273|
| YOR370C   | No stall |  -1.273|
| YAL031C   | Stall    |  -1.275|
| YHL010C   | Stall    |  -1.277|
| YAL058W   | No stall |  -1.278|
| YGL133W   | Stall    |  -1.278|
| YDR186C   | No stall |  -1.280|
| YDR137W   | No stall |  -1.282|
| YOR140W   | No stall |  -1.282|
| YGL167C   | Stall    |  -1.283|
| YKR028W   | No stall |  -1.284|
| YOL082W   | No stall |  -1.284|
| YPR122W   | Stall    |  -1.285|
| YPL208W   | No stall |  -1.285|
| YDR108W   | No stall |  -1.287|
| YNL056W   | No stall |  -1.288|
| YPL082C   | Stall    |  -1.288|
| YOL013C   | Stall    |  -1.288|
| YOL130W   | No stall |  -1.289|
| YNL025C   | No stall |  -1.290|
| YGL126W   | No stall |  -1.292|
| YOR355W   | No stall |  -1.294|
| YHR073W   | Stall    |  -1.294|
| YLR247C   | Stall    |  -1.297|
| YIL126W   | Stall    |  -1.301|
| YOR354C   | No stall |  -1.302|
| YLR099C   | No stall |  -1.303|
| YKL219W   | No stall |  -1.304|
| YOR213C   | No stall |  -1.304|
| YHR077C   | No stall |  -1.304|
| YNR055C   | No stall |  -1.306|
| YNL238W   | Stall    |  -1.309|
| YDR028C   | Stall    |  -1.311|
| YPL011C   | No stall |  -1.312|
| YKL043W   | No stall |  -1.316|
| YMR162C   | No stall |  -1.316|
| YOR002W   | No stall |  -1.317|
| YEL052W   | No stall |  -1.317|
| YDR180W   | Stall    |  -1.323|
| YBL008W   | No stall |  -1.323|
| YEL006W   | No stall |  -1.324|
| YJR046W   | No stall |  -1.329|
| YKL203C   | No stall |  -1.329|
| YMR191W   | No stall |  -1.330|
| YML117W   | Stall    |  -1.331|
| YGL006W   | Stall    |  -1.331|
| YOR324C   | No stall |  -1.338|
| YPL115C   | Stall    |  -1.338|
| YNL236W   | No stall |  -1.339|
| YJL076W   | Stall    |  -1.340|
| YDR325W   | Stall    |  -1.340|
| YGR223C   | No stall |  -1.348|
| YDL013W   | Stall    |  -1.350|
| YNL186W   | Stall    |  -1.351|
| YMR026C   | No stall |  -1.351|
| YOR118W   | No stall |  -1.351|
| YPL134C   | No stall |  -1.351|
| YPL193W   | Stall    |  -1.351|
| YDR288W   | No stall |  -1.352|
| YNL261W   | Stall    |  -1.356|
| YPL207W   | Stall    |  -1.358|
| YHR103W   | No stall |  -1.358|
| YKR099W   | Stall    |  -1.360|
| YKL101W   | No stall |  -1.361|
| YLR115W   | Stall    |  -1.362|
| YPL122C   | Stall    |  -1.363|
| YJR137C   | Stall    |  -1.364|
| YLR168C   | No stall |  -1.369|
| YNL053W   | No stall |  -1.373|
| YML046W   | No stall |  -1.374|
| YHR120W   | No stall |  -1.375|
| YOR222W   | No stall |  -1.378|
| YER169W   | Stall    |  -1.382|
| YJR122W   | No stall |  -1.383|
| YNL023C   | No stall |  -1.385|
| YLR427W   | Stall    |  -1.386|
| YPR067W   | No stall |  -1.388|
| YBR153W   | No stall |  -1.391|
| YKL197C   | No stall |  -1.392|
| YBL105C   | Stall    |  -1.393|
| YOR017W   | No stall |  -1.393|
| YJL109C   | Stall    |  -1.393|
| YGL233W   | No stall |  -1.397|
| YCL010C   | No stall |  -1.399|
| YBR133C   | No stall |  -1.399|
| YER111C   | Stall    |  -1.401|
| YDR389W   | Stall    |  -1.402|
| YMR109W   | Stall    |  -1.403|
| YBR183W   | No stall |  -1.405|
| YKR010C   | Stall    |  -1.407|
| YER129W   | No stall |  -1.409|
| YDR217C   | Stall    |  -1.409|
| YMR020W   | No stall |  -1.415|
| YPL137C   | No stall |  -1.417|
| YLR073C   | No stall |  -1.421|
| YAL014C   | No stall |  -1.421|
| YDR138W   | Stall    |  -1.422|
| YOR329C   | Stall    |  -1.423|
| YIL049W   | No stall |  -1.428|
| YLR023C   | No stall |  -1.430|
| YAR003W   | No stall |  -1.434|
| YGR261C   | Stall    |  -1.434|
| YKL098W   | No stall |  -1.436|
| YPR091C   | Stall    |  -1.437|
| YOL136C   | No stall |  -1.440|
| YDR313C   | No stall |  -1.441|
| YJL059W   | No stall |  -1.443|
| YCL051W   | Stall    |  -1.446|
| YOR191W   | Stall    |  -1.446|
| YOR043W   | Stall    |  -1.449|
| YKL072W   | Stall    |  -1.450|
| YMR313C   | No stall |  -1.456|
| YDR407C   | No stall |  -1.456|
| YHR067W   | No stall |  -1.459|
| YLR383W   | No stall |  -1.461|
| YDR176W   | Stall    |  -1.461|
| YOR326W   | Stall    |  -1.463|
| YOR005C   | Stall    |  -1.466|
| YMR078C   | Stall    |  -1.467|
| YGR157W   | No stall |  -1.467|
| YPL003W   | No stall |  -1.467|
| YML038C   | No stall |  -1.469|
| YDR251W   | Stall    |  -1.471|
| YCL005W   | No stall |  -1.473|
| YOR303W   | No stall |  -1.476|
| YHL023C   | Stall    |  -1.476|
| YLR086W   | Stall    |  -1.483|
| YPR105C   | No stall |  -1.488|
| YDR488C   | No stall |  -1.491|
| YLR114C   | No stall |  -1.491|
| YAL009W   | No stall |  -1.491|
| YOL069W   | No stall |  -1.492|
| YKR022C   | No stall |  -1.493|
| YFL027C   | No stall |  -1.495|
| YKL194C   | No stall |  -1.497|
| YER152C   | No stall |  -1.497|
| YGR024C   | No stall |  -1.498|
| YGL155W   | No stall |  -1.501|
| YAL040C   | No stall |  -1.505|
| YHR129C   | No stall |  -1.505|
| YJL095W   | No stall |  -1.505|
| YDR038C   | No stall |  -1.513|
| YDR039C   | No stall |  -1.513|
| YBR158W   | No stall |  -1.515|
| YDR466W   | Stall    |  -1.516|
| YLR007W   | No stall |  -1.519|
| YKR096W   | Stall    |  -1.519|
| YJL044C   | No stall |  -1.519|
| YPL221W   | No stall |  -1.519|
| YOR127W   | Stall    |  -1.522|
| YPL214C   | No stall |  -1.524|
| YDR011W   | Stall    |  -1.525|
| YIL002C   | No stall |  -1.527|
| YPR164W   | No stall |  -1.530|
| YKL033W   | No stall |  -1.532|
| YMR171C   | No stall |  -1.532|
| YOL062C   | No stall |  -1.538|
| YMR127C   | No stall |  -1.540|
| YHR106W   | No stall |  -1.543|
| YMR179W   | No stall |  -1.546|
| YJL083W   | Stall    |  -1.546|
| YGL236C   | No stall |  -1.548|
| YKR021W   | Stall    |  -1.550|
| YLR188W   | No stall |  -1.554|
| YJL046W   | No stall |  -1.554|
| YDR334W   | Stall    |  -1.558|
| YER144C   | No stall |  -1.561|
| YGR258C   | Stall    |  -1.563|
| YIR002C   | Stall    |  -1.565|
| YGL247W   | No stall |  -1.572|
| YPR138C   | No stall |  -1.573|
| YNL188W   | No stall |  -1.574|
| YDR524C   | No stall |  -1.575|
| YLR272C   | No stall |  -1.575|
| YIL159W   | Stall    |  -1.582|
| YCR093W   | No stall |  -1.583|
| YBR030W   | No stall |  -1.584|
| YLR265C   | Stall    |  -1.588|
| YOR026W   | No stall |  -1.592|
| YJL197W   | Stall    |  -1.592|
| YJL054W   | No stall |  -1.593|
| YLR084C   | Stall    |  -1.598|
| YDR369C   | Stall    |  -1.601|
| YMR218C   | Stall    |  -1.602|
| YDR349C   | No stall |  -1.607|
| YLR375W   | Stall    |  -1.607|
| YHL047C   | No stall |  -1.608|
| YLR353W   | No stall |  -1.609|
| YJR035W   | Stall    |  -1.609|
| YIL149C   | No stall |  -1.613|
| YOR231W   | No stall |  -1.614|
| YLR360W   | No stall |  -1.615|
| YGR077C   | No stall |  -1.615|
| YDR393W   | No stall |  -1.620|
| YOL100W   | Stall    |  -1.620|
| YIL156W   | Stall    |  -1.622|
| YNL101W   | No stall |  -1.622|
| YMR066W   | Stall    |  -1.624|
| YBR140C   | No stall |  -1.625|
| YLR412W   | Stall    |  -1.626|
| YKR079C   | No stall |  -1.628|
| YKL092C   | Stall    |  -1.633|
| YGL094C   | No stall |  -1.638|
| YJL071W   | Stall    |  -1.646|
| YCL024W   | Stall    |  -1.648|
| YNL273W   | Stall    |  -1.649|
| YPL024W   | No stall |  -1.654|
| YBR274W   | No stall |  -1.655|
| YPL015C   | No stall |  -1.655|
| YGR222W   | No stall |  -1.657|
| YGR089W   | Stall    |  -1.658|
| YCR033W   | Stall    |  -1.658|
| YMR172W   | Stall    |  -1.659|
| YNL283C   | No stall |  -1.662|
| YPL179W   | Stall    |  -1.663|
| YER005W   | Stall    |  -1.664|
| YDR310C   | Stall    |  -1.664|
| YBR228W   | No stall |  -1.668|
| YPL064C   | No stall |  -1.671|
| YPL040C   | No stall |  -1.671|
| YPR009W   | Stall    |  -1.673|
| YOR346W   | No stall |  -1.676|
| YLL019C   | Stall    |  -1.687|
| YNL063W   | No stall |  -1.692|
| YLR289W   | Stall    |  -1.692|
| YOR336W   | No stall |  -1.695|
| YEL043W   | No stall |  -1.699|
| YNL073W   | No stall |  -1.700|
| YOR306C   | No stall |  -1.702|
| YHR119W   | Stall    |  -1.704|
| YJR138W   | Stall    |  -1.706|
| YOR275C   | No stall |  -1.711|
| YAL056W   | No stall |  -1.711|
| YOL140W   | No stall |  -1.714|
| YER060W   | No stall |  -1.718|
| YOR073W   | Stall    |  -1.721|
| YPL006W   | No stall |  -1.721|
| YLR097C   | No stall |  -1.725|
| YBL011W   | Stall    |  -1.739|
| YDR528W   | Stall    |  -1.740|
| YOL063C   | No stall |  -1.741|
| YER081W   | No stall |  -1.742|
| YOL081W   | No stall |  -1.748|
| YDR443C   | No stall |  -1.755|
| YOR149C   | No stall |  -1.757|
| YOR034C   | No stall |  -1.762|
| YNL329C   | No stall |  -1.765|
| YKR031C   | Stall    |  -1.767|
| YGR184C   | Stall    |  -1.772|
| YGL257C   | No stall |  -1.772|
| YPR095C   | Stall    |  -1.778|
| YOR334W   | Stall    |  -1.783|
| YOR064C   | Stall    |  -1.783|
| YOR291W   | Stall    |  -1.784|
| YDR359C   | No stall |  -1.793|
| YFL007W   | Stall    |  -1.794|
| YDR073W   | No stall |  -1.794|
| YGR057C   | No stall |  -1.795|
| YLR443W   | No stall |  -1.799|
| YNL257C   | No stall |  -1.800|
| YNL272C   | Stall    |  -1.801|
| YMR282C   | No stall |  -1.805|
| YMR100W   | Stall    |  -1.809|
| YGL197W   | Stall    |  -1.811|
| YLL051C   | No stall |  -1.816|
| YNL038W   | No stall |  -1.819|
| YIL085C   | No stall |  -1.821|
| YMR176W   | No stall |  -1.829|
| YPR057W   | No stall |  -1.830|
| YJL210W   | No stall |  -1.833|
| YGL136C   | No stall |  -1.833|
| YMR119W   | No stall |  -1.838|
| YKL048C   | Stall    |  -1.838|
| YKR054C   | No stall |  -1.845|
| YKR044W   | Stall    |  -1.872|
| YNL080C   | Stall    |  -1.873|
| YKL015W   | No stall |  -1.873|
| YCR039C   | No stall |  -1.873|
| YJL187C   | Stall    |  -1.879|
| YDR040C   | No stall |  -1.885|
| YDR128W   | No stall |  -1.885|
| YBR275C   | Stall    |  -1.893|
| YKR036C   | No stall |  -1.896|
| YDR326C   | Stall    |  -1.896|
| YCR008W   | No stall |  -1.901|
| YMR299C   | No stall |  -1.907|
| YOR330C   | Stall    |  -1.907|
| YPR115W   | No stall |  -1.915|
| YHR204W   | No stall |  -1.918|
| YGR188C   | Stall    |  -1.942|
| YIL139C   | No stall |  -1.943|
| YDR283C   | Stall    |  -1.943|
| YNL291C   | No stall |  -1.948|
| YLR138W   | Stall    |  -1.949|
| YFL033C   | Stall    |  -1.950|
| YLR189C   | Stall    |  -1.955|
| YPL244C   | No stall |  -1.969|
| YDL146W   | Stall    |  -1.973|
| YLR310C   | Stall    |  -1.979|
| YDL080C   | No stall |  -1.982|
| YER098W   | Stall    |  -1.989|
| YLR305C   | Stall    |  -1.994|
| YLR130C   | No stall |  -1.997|
| YJR089W   | Stall    |  -1.999|
| YGL113W   | Stall    |  -2.000|
| YLR368W   | No stall |  -2.002|
| YBR168W   | Stall    |  -2.005|
| YLR047C   | No stall |  -2.021|
| YGL134W   | No stall |  -2.030|
| YHR155W   | No stall |  -2.031|
| YKR086W   | Stall    |  -2.064|
| YMR270C   | Stall    |  -2.069|
| YDR206W   | No stall |  -2.079|
| YLR430W   | Stall    |  -2.083|
| YNL029C   | Stall    |  -2.087|
| YAL029C   | Stall    |  -2.098|
| YIL147C   | No stall |  -2.101|
| YOR080W   | Stall    |  -2.108|
| YMR180C   | No stall |  -2.109|
| YHR082C   | Stall    |  -2.113|
| YIL031W   | No stall |  -2.119|
| YOL103W   | No stall |  -2.123|
| YJL093C   | No stall |  -2.125|
| YEL016C   | No stall |  -2.138|
| YBL004W   | Stall    |  -2.139|
| YNL271C   | Stall    |  -2.146|
| YOR205C   | No stall |  -2.149|
| YDL238C   | No stall |  -2.152|
| YER155C   | Stall    |  -2.162|
| YPL005W   | Stall    |  -2.174|
| YCR023C   | No stall |  -2.188|
| YOR305W   | No stall |  -2.189|
| YIL146C   | Stall    |  -2.196|
| YHR099W   | Stall    |  -2.202|
| YDR285W   | Stall    |  -2.228|
| YPL149W   | No stall |  -2.244|
| YKL201C   | Stall    |  -2.263|
| YLR024C   | Stall    |  -2.269|
| YER153C   | No stall |  -2.305|
| YLL049W   | No stall |  -2.312|
| YPR021C   | Stall    |  -2.322|
| YGR257C   | No stall |  -2.323|
| YGL256W   | No stall |  -2.341|
| YMR064W   | No stall |  -2.356|
| YHR079C   | Stall    |  -2.360|
| YML023C   | No stall |  -2.360|
| YER173W   | Stall    |  -2.365|
| YKR098C   | Stall    |  -2.373|
| YOR193W   | No stall |  -2.377|
| YIL173W   | Stall    |  -2.382|
| YJL222W   | Stall    |  -2.382|
| YGL215W   | No stall |  -2.433|
| YDR545W   | No stall |  -2.462|
| YLR467W   | No stall |  -2.462|
| YJL209W   | No stall |  -2.469|
| YPL022W   | Stall    |  -2.472|
| YHR151C   | No stall |  -2.495|
| YJL092W   | Stall    |  -2.495|
| YPL072W   | Stall    |  -2.497|
| YDR420W   | Stall    |  -2.512|
| YBR163W   | No stall |  -2.543|
| YKL176C   | No stall |  -2.544|
| YGR296W   | No stall |  -2.554|
| YGL229C   | Stall    |  -2.561|
| YMR207C   | No stall |  -2.568|
| YOL068C   | Stall    |  -2.584|
| YGL131C   | Stall    |  -2.604|
| YOR030W   | Stall    |  -2.623|
| YNL294C   | Stall    |  -2.624|
| YGR098C   | Stall    |  -2.646|
| YNL339C   | Stall    |  -2.679|
| YPL283C   | Stall    |  -2.679|
| YMR190C   | Stall    |  -2.763|
| YBL014C   | Stall    |  -2.787|
| YLR320W   | Stall    |  -2.799|
| YOL011W   | No stall |  -2.802|
| YDR459C   | Stall    |  -2.838|
| YNL094W   | Stall    |  -2.841|
| YJR147W   | No stall |  -2.880|
| YPL075W   | No stall |  -2.890|
| YDL138W   | Stall    |  -2.897|
| YML104C   | No stall |  -2.901|
| YDL197C   | No stall |  -2.939|
| YBL052C   | Stall    |  -2.956|
| YER132C   | Stall    |  -2.961|
| YBR276C   | No stall |  -3.014|
| YER041W   | Stall    |  -3.021|
| YLR466W   | No stall |  -3.184|
| YER190W   | No stall |  -3.197|
| YFL060C   | No stall |  -3.223|
| YHR071W   | No stall |  -3.235|
| YNL334C   | No stall |  -3.236|
| YJL129C   | Stall    |  -3.260|
| YPL174C   | Stall    |  -3.279|
| YJR120W   | No stall |  -3.297|
| YHL032C   | No stall |  -3.303|
| YLR087C   | Stall    |  -3.306|
| YOL113W   | No stall |  -3.465|
| YDR275W   | No stall |  -3.662|
| YHL035C   | No stall |  -3.726|
| YML051W   | No stall |  -4.012|
| YEL009C   | Stall    |  -5.320|
