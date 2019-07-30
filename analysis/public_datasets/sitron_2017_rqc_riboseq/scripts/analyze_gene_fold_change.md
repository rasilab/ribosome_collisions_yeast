Analyze changes in gene ribosome density between WT and HEL2/ASC1 mutants
================
rasi
30 July, 2019

-   [Import libraries](#import-libraries)
-   [Load SGD annotations](#load-sgd-annotations)
-   [Get transcript annotations](#get-transcript-annotations)
-   [Load RQC stalls](#load-rqc-stalls)
-   [Read count data](#read-count-data)
-   [Check that the genotypes are correct](#check-that-the-genotypes-are-correct)
-   [Sum only read counts for gene that have &gt; `count_threshold` reads in all samples](#sum-only-read-counts-for-gene-that-have-count_threshold-reads-in-all-samples)
-   [Prepare column data for DESeq 2 input](#prepare-column-data-for-deseq-2-input)
-   [Run DESeq2](#run-deseq2)
-   [Calculate log2 fold-changes between WT and mutant strains](#calculate-log2-fold-changes-between-wt-and-mutant-strains)
-   [Look at genes that are up-regulated in asc1 and hel2 KO](#look-at-genes-that-are-up-regulated-in-asc1-and-hel2-ko)
-   [Join fold-change with RQC stall presence](#join-fold-change-with-rqc-stall-presence)
-   [Plot log2 fold-change ASC1 / HEL2 vs WT as a function of stall strength](#plot-log2-fold-change-asc1-hel2-vs-wt-as-a-function-of-stall-strength)
-   [Test if stall-containing genes have lower or higher Log2 fold change between ASC1+HEL2 KO vs WT](#test-if-stall-containing-genes-have-lower-or-higher-log2-fold-change-between-asc1hel2-ko-vs-wt)
-   [Source data for Fig 6D](#source-data-for-fig-6d)

Import libraries
================

``` r
library(data.table)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(org.Sc.sgd.db)
library(biobroom)
library(DESeq2)
library(tidyverse)
library(rasilabRtemplates)

# genes having below these counts in any sample are discarded
count_threshold <- 100
```

Load SGD annotations
====================

``` r
gene_annotations <- "/fh/fast/subramaniam_a/db/rasi/genomes/yeast/Saccharomyces_cerevisiae/sgd/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff" %>% 
  rtracklayer::readGFF() %>% 
  as_tibble() %>% 
  filter(type == "gene" & orf_classification == "Verified") %>% 
  select(gene, ID, Note) %>%
  print()
```

    ## # A tibble: 4,896 x 3
    ##    gene  ID      Note     
    ##    <chr> <chr>   <I(list)>
    ##  1 PAU8  YAL068C <chr [1]>
    ##  2 SEO1  YAL067C <chr [1]>
    ##  3 <NA>  YAL064W <chr [1]>
    ##  4 FLO9  YAL063C <chr [1]>
    ##  5 GDH3  YAL062W <chr [1]>
    ##  6 BDH1  YAL060W <chr [1]>
    ##  7 ECM1  YAL059W <chr [1]>
    ##  8 CNE1  YAL058W <chr [1]>
    ##  9 GPB2  YAL056W <chr [1]>
    ## 10 PEX22 YAL055W <chr [1]>
    ## # ... with 4,886 more rows

Get transcript annotations
==========================

``` r
tx <- transcripts(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

tx_annotations <- org.Sc.sgd.db %>% 
  AnnotationDbi::select(keys = keys(., keytype = "ENSEMBLTRANS"), 
         keytype = "ENSEMBLTRANS", 
         columns = c("ENSEMBLTRANS", "GENENAME", "DESCRIPTION")) %>% 
  as_tibble() %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  right_join(tidy(tx) %>% dplyr::select(tx_id, tx_name, strand), 
             by = c("ensembltrans" = "tx_name")) %>% 
  group_by(tx_id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  print()
```

    ## # A tibble: 6,692 x 5
    ##    ensembltrans genename description                          tx_id strand
    ##    <chr>        <chr>    <chr>                                <int> <chr> 
    ##  1 YAL069W      <NA>     <NA>                                     1 +     
    ##  2 YAL068W-A    <NA>     <NA>                                     2 +     
    ##  3 YAL067W-A    <NA>     Putative protein of unknown functio…     3 +     
    ##  4 YAL066W      <NA>     <NA>                                     4 +     
    ##  5 YAL064W-B    <NA>     Fungal-specific protein of unknown …     5 +     
    ##  6 YAL064W      <NA>     Protein of unknown function; may in…     6 +     
    ##  7 YAL062W      GDH3     NADP(+)-dependent glutamate dehydro…     7 +     
    ##  8 YAL061W      BDH2     Putative medium-chain alcohol dehyd…     8 +     
    ##  9 YAL060W      BDH1     NAD-dependent (R,R)-butanediol dehy…     9 +     
    ## 10 YAL059W      ECM1     Pre-ribosomal factor involved in 60…    10 +     
    ## # ... with 6,682 more rows

Load RQC stalls
===============

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

Read count data
===============

``` r
counts <- list.files("../processeddata/", pattern = "tx_read_counts.tsv", 
                    recursive = T, full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(sample = str_extract(file, "[^/]+(?=/+tx_read_counts.tsv)")) %>% 
  separate(sample, c("genotype", "reporter", "sampletype"), remove = F) %>% 
  mutate(data = map(file, read_tsv)) %>%
  select(-sno, -file) %>% 
  unnest() %>% 
  print()
```

    ## # A tibble: 25,964 x 7
    ##    sample             genotype reporter sampletype txHit strand alncount
    ##    <chr>              <chr>    <chr>    <chr>      <int> <chr>     <int>
    ##  1 asc1_nonstall_mrna asc1     nonstall mrna          NA -       2043844
    ##  2 asc1_nonstall_mrna asc1     nonstall mrna          NA +       1084530
    ##  3 asc1_nonstall_mrna asc1     nonstall mrna        4146 +         79313
    ##  4 asc1_nonstall_mrna asc1     nonstall mrna        6337 +         71829
    ##  5 asc1_nonstall_mrna asc1     nonstall mrna        5576 +         71317
    ##  6 asc1_nonstall_mrna asc1     nonstall mrna         235 +         70537
    ##  7 asc1_nonstall_mrna asc1     nonstall mrna        4643 +         70498
    ##  8 asc1_nonstall_mrna asc1     nonstall mrna        3383 +         69837
    ##  9 asc1_nonstall_mrna asc1     nonstall mrna        4640 +         69526
    ## 10 asc1_nonstall_mrna asc1     nonstall mrna        3385 +         68613
    ## # ... with 25,954 more rows

Check that the genotypes are correct
====================================

We see that the genes that are deleted have very low counts as expected.

``` r
genotype_genes <- c(
  "HEL2" = "YDR266C",
  "ASC1" = "YMR116C",
  "SLH1" = "YGR271W")

counts %>% 
  left_join(tx_annotations, by = c("txHit" = "tx_id")) %>%
  filter(ensembltrans %in% genotype_genes) %>%
  arrange(desc(genename)) %>% 
  select(genename, sample, genotype, alncount) %>%
  knitr::kable()
```

| genename | sample               | genotype |  alncount|
|:---------|:---------------------|:---------|---------:|
| SLH1     | asc1\_nonstall\_mrna | asc1     |      1353|
| SLH1     | hel2\_nonstall\_mrna | hel2     |      1214|
| SLH1     | slh1\_nonstall\_mrna | slh1     |         7|
| SLH1     | wt\_nonstall\_mrna   | wt       |      1268|
| HEL2     | asc1\_nonstall\_mrna | asc1     |      1096|
| HEL2     | hel2\_nonstall\_mrna | hel2     |         9|
| HEL2     | slh1\_nonstall\_mrna | slh1     |      1110|
| HEL2     | wt\_nonstall\_mrna   | wt       |       799|
| ASC1     | asc1\_nonstall\_mrna | asc1     |        32|
| ASC1     | hel2\_nonstall\_mrna | hel2     |      4174|
| ASC1     | slh1\_nonstall\_mrna | slh1     |      4951|
| ASC1     | wt\_nonstall\_mrna   | wt       |      3728|

Sum only read counts for gene that have &gt; `count_threshold` reads in all samples
===================================================================================

``` r
count_data <- counts %>% 
  group_by(txHit, sample) %>%
  summarize(count = sum(alncount)) %>%
  ungroup() %>%
  filter(!is.na(txHit)) %>%
  filter(count > count_threshold) %>% 
  filter(str_detect(sample, "nonstall")) %>%
  filter(str_detect(sample, "mrna")) %>%
  select(sample, txHit, count) %>%
  spread(sample, count) %>%
  filter_all(all_vars(!is.na(.))) %>%
  inner_join(tx_annotations %>% select(tx_id, ensembltrans), by = c("txHit" = "tx_id")) %>%
  inner_join(gene_annotations %>% select(ID), by = c("ensembltrans" = "ID")) %>%
  select(-ensembltrans) %>% 
  as.data.frame() %>%
  column_to_rownames(var = "txHit")

as_tibble(count_data)
```

    ## # A tibble: 4,354 x 4
    ##    asc1_nonstall_mrna hel2_nonstall_mr… slh1_nonstall_mr… wt_nonstall_mrna
    ##  *              <int>             <int>             <int>            <int>
    ##  1                723               546               455              649
    ##  2                704              1282              2028              922
    ##  3                630               453               594              431
    ##  4               1205              1024              1292             1304
    ##  5                247               131               156              228
    ##  6               1290              2040              2134             1067
    ##  7                908               849               885              768
    ##  8               4027              2958              3796             3418
    ##  9               1406              1197              1414             1160
    ## 10              48431             37266             41673            29905
    ## # ... with 4,344 more rows

Prepare column data for DESeq 2 input
=====================================

``` r
col_data <- colnames(count_data) %>% 
  enframe("sno", "sample") %>% 
  separate(sample, c("genotype", "reporter", "sampletype"), remove = F) %>% 
  mutate(genotype = if_else(genotype %in% c("asc1", "hel2"), "asc1hel2", genotype)) %>%
  select(-sno) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "sample") %>% 
  # set WT to be the reference level for log2 fold-changes in DESeq2
  mutate(genotype = relevel(factor(genotype), ref = "wt"))

col_data
```

    ##   genotype reporter sampletype
    ## 1 asc1hel2 nonstall       mrna
    ## 2 asc1hel2 nonstall       mrna
    ## 3     slh1 nonstall       mrna
    ## 4       wt nonstall       mrna

Run DESeq2
==========

``` r
dds <- DESeqDataSetFromMatrix(count_data, col_data, ~ genotype)
dds <- DESeq(dds)
```

Calculate log2 fold-changes between WT and mutant strains
=========================================================

``` r
# skip the Intercept column
lfc <- resultsNames(dds)[2:3] %>% 
  enframe() %>% 
  mutate(deseq_results = map(value, function(x) DESeq2::results(dds, contrast = list(x)))) %>%
  mutate(lfc = map(deseq_results, function(res)
    res[c('log2FoldChange', 'baseMean')] %>% as.data.frame() %>% rownames_to_column("tx_id"))) %>%
  mutate(samplepair = str_extract(value, "(?<=genotype_).+")) %>%
  select(-deseq_results, -name, -value) %>%
  unnest() %>%
  rename(lfc = log2FoldChange) %>%
  spread(samplepair, lfc) %>%
  mutate(tx_id = as.integer(tx_id)) %>% 
  left_join(tx_annotations, by = c("tx_id")) %>%
  rename(id = ensembltrans) %>% 
  print()  
```

    ## # A tibble: 4,354 x 8
    ##    tx_id baseMean asc1hel2_vs_wt slh1_vs_wt id    genename description
    ##    <int>    <dbl>          <dbl>      <dbl> <chr> <chr>    <chr>      
    ##  1    10    1245.         0.127      0.958  YAL0… ECM1     Pre-riboso…
    ##  2   100     469.         0.446      0.0724 YAL0… SYN8     Endosomal …
    ##  3  1000     493.        -0.445     -0.638  YDR2… COQ4     Protein wi…
    ##  4  1001     503.         0.428      0.449  YDR2… MSC2     Endoplasmi…
    ##  5  1002    1556.         0.437      0.636  YDR2… EBS1     Protein in…
    ##  6  1003    1592.         0.255      0.419  YDR2… MSS4     Phosphatid…
    ##  7  1007    3415.        -0.0275     0.206  YDR2… GCD6     Catalytic …
    ##  8  1008    2877.        -0.252     -0.237  YDR2… TCP1     Alpha subu…
    ##  9  1009     794.        -0.0714     0.185  YDR2… UPC2     Sterol reg…
    ## 10   101     334.         0.354      0.558  YAL0… MDM10    Subunit of…
    ## # ... with 4,344 more rows, and 1 more variable: strand <chr>

Look at genes that are up-regulated in asc1 and hel2 KO
=======================================================

``` r
lfc %>% 
  arrange(desc(asc1hel2_vs_wt)) %>%
  select(asc1hel2_vs_wt, genename, everything(), -tx_id, description, -strand) %>%
  print()
```

    ## # A tibble: 4,354 x 6
    ##    asc1hel2_vs_wt genename baseMean slh1_vs_wt id     description         
    ##             <dbl> <chr>       <dbl>      <dbl> <chr>  <chr>               
    ##  1           3.36 PRM7        1745.      3.62  YDL03… Pheromone-regulated…
    ##  2           2.96 FLR1        1187.      2.17  YBR00… Plasma membrane tra…
    ##  3           2.68 AGA1         709.      2.37  YNR04… Anchorage subunit o…
    ##  4           2.51 GRE2        1919.      2.03  YOL15… 3-methylbutanal red…
    ##  5           2.31 CHA1        1251.      2.58  YCL06… Catabolic L-serine …
    ##  6           2.19 HMS2         650.      1.75  YJR14… Protein with simila…
    ##  7           2.18 HSP12       1218.      0.155 YFL01… Plasma membrane pro…
    ##  8           2.14 SPI1         496.      0.596 YER15… GPI-anchored cell w…
    ##  9           2.03 <NA>         681.      1.68  YOR30… CPA1 uORF; Arginine…
    ## 10           1.91 BUD31        356.      1.52  YCR06… Component of the SF…
    ## # ... with 4,344 more rows

Join fold-change with RQC stall presence
========================================

``` r
lfc_stall_data <- lfc %>% 
  left_join(rqc_stalls, by = "id") %>% 
  arrange(desc(asc1hel2_vs_wt)) %>% 
  print()
```

    ## # A tibble: 4,354 x 11
    ##    tx_id baseMean asc1hel2_vs_wt slh1_vs_wt id    genename description
    ##    <int>    <dbl>          <dbl>      <dbl> <chr> <chr>    <chr>      
    ##  1  1295    1745.           3.36      3.62  YDL0… PRM7     Pheromone-…
    ##  2   406    1187.           2.96      2.17  YBR0… FLR1     Plasma mem…
    ##  3  5321     709.           2.68      2.37  YNR0… AGA1     Anchorage …
    ##  4  5550    1919.           2.51      2.03  YOL1… GRE2     3-methylbu…
    ##  5   667    1251.           2.31      2.58  YCL0… CHA1     Catabolic …
    ##  6  3452     650.           2.19      1.75  YJR1… HMS2     Protein wi…
    ##  7  1975    1218.           2.18      0.155 YFL0… HSP12    Plasma mem…
    ##  8  1758     496.           2.14      0.596 YER1… SPI1     GPI-anchor…
    ##  9  5803     681.           2.03      1.68  YOR3… <NA>     CPA1 uORF;…
    ## 10   637     356.           1.91      1.52  YCR0… BUD31    Component …
    ## # ... with 4,344 more rows, and 4 more variables: strand <chr>, pos <int>,
    ## #   ngram_weight <int>, ngram <chr>

Plot log2 fold-change ASC1 / HEL2 vs WT as a function of stall strength
=======================================================================

``` r
plot_data <- lfc_stall_data %>% 
  mutate(ngram_weight = as.factor(if_else(is.na(ngram_weight), 0, 1))) %>% 
  group_by(ngram_weight) %>% 
  mutate(`n` = paste0("N = ", dplyr::n())) %>% 
  ungroup() %>% 
  mutate(ngram_weight = fct_recode(ngram_weight, `No stall` = "0", `Stall` = "1"))
  
plot_data %>% 
  ggplot(aes(x = ngram_weight, y = asc1hel2_vs_wt, fill = ngram_weight)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  labs(x = "S. cerevisiae genes", y = "ASC1 & HEL2 / WT (log2, a.u.)") +
  geom_text(aes(x = ngram_weight, label = n),
            data = plot_data %>% group_by(ngram_weight) %>% slice(1),
            y = -2.1, size = 2.8) +
  scale_y_continuous(limits = c(-2.1, NA)) +
  scale_fill_manual(values = cbPalette, guide = 'none') + 
  NULL
```

![](analyze_gene_fold_change_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
ggsave("../figures/distribution_of_asc1hel2ko_lfc_for_rqc_stall_containing_saccer_genes.pdf")
```

Test if stall-containing genes have lower or higher Log2 fold change between ASC1+HEL2 KO vs WT
===============================================================================================

``` r
wilcox.test(asc1hel2_vs_wt ~ ngram_weight, data = plot_data, alternative = "two.sided")
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  asc1hel2_vs_wt by ngram_weight
    ## W = 1714200, p-value = 0.7003
    ## alternative hypothesis: true location shift is not equal to 0

Source data for Fig 6D
=======================

``` r
plot_data %>% 
  select(id, genename, ngram_weight, asc1hel2_vs_wt) %>% 
  rename(x = ngram_weight, y = asc1hel2_vs_wt) %>% 
  mutate_if(is.numeric, funs(round(., 3))) %>% 
  knitr::kable()
```

| id        | genename  | x        |       y|
|:----------|:----------|:---------|-------:|
| YDL039C   | PRM7      | No stall |   3.356|
| YBR008C   | FLR1      | No stall |   2.959|
| YNR044W   | AGA1      | No stall |   2.679|
| YOL151W   | GRE2      | No stall |   2.510|
| YCL064C   | CHA1      | No stall |   2.313|
| YJR147W   | HMS2      | No stall |   2.187|
| YFL014W   | HSP12     | No stall |   2.179|
| YER150W   | SPI1      | No stall |   2.142|
| YOR302W   | NA        | No stall |   2.028|
| YCR063W   | BUD31     | No stall |   1.909|
| YHR096C   | HXT5      | No stall |   1.906|
| YEL021W   | URA3      | No stall |   1.868|
| YGR088W   | CTT1      | Stall    |   1.760|
| YOR306C   | MCH5      | No stall |   1.757|
| YKR013W   | PRY2      | No stall |   1.732|
| YCL037C   | SRO9      | No stall |   1.724|
| YBR010W   | HHT2      | No stall |   1.697|
| YCR089W   | FIG2      | No stall |   1.664|
| YBR072W   | HSP26     | No stall |   1.661|
| YGR268C   | HUA1      | No stall |   1.651|
| YFL056C   | AAD6      | No stall |   1.642|
| YJL078C   | PRY3      | No stall |   1.636|
| YDL123W   | SNA4      | No stall |   1.636|
| YOR247W   | SRL1      | No stall |   1.595|
| YNL160W   | YGP1      | No stall |   1.534|
| YGR279C   | SCW4      | No stall |   1.527|
| YHR070W   | TRM5      | No stall |   1.510|
| YER103W   | SSA4      | No stall |   1.504|
| YLL018C-A | COX19     | No stall |   1.489|
| YGL076C   | RPL7A     | No stall |   1.473|
| YLR337C   | VRP1      | Stall    |   1.470|
| YBR016W   | NA        | No stall |   1.458|
| YHR143W   | DSE2      | No stall |   1.428|
| YML116W   | ATR1      | No stall |   1.418|
| YPR154W   | PIN3      | No stall |   1.392|
| YDL204W   | RTN2      | No stall |   1.370|
| YJR070C   | LIA1      | No stall |   1.349|
| YDR126W   | SWF1      | No stall |   1.346|
| YJR077C   | MIR1      | No stall |   1.331|
| YAL034C   | FUN19     | No stall |   1.327|
| YMR096W   | SNZ1      | No stall |   1.323|
| YHR025W   | THR1      | No stall |   1.293|
| YER045C   | ACA1      | No stall |   1.288|
| YOL105C   | WSC3      | No stall |   1.278|
| YBR135W   | CKS1      | No stall |   1.258|
| YPR183W   | DPM1      | No stall |   1.247|
| YGL128C   | CWC23     | No stall |   1.245|
| YMR169C   | ALD3      | No stall |   1.236|
| YLR073C   | RFU1      | No stall |   1.226|
| YKL018W   | SWD2      | No stall |   1.226|
| YJL079C   | PRY1      | No stall |   1.225|
| YNL096C   | RPS7B     | No stall |   1.216|
| YDR432W   | NPL3      | Stall    |   1.207|
| YLR263W   | RED1      | Stall    |   1.202|
| YCL005W   | LDB16     | No stall |   1.195|
| YLR286C   | CTS1      | No stall |   1.185|
| YLR438W   | CAR2      | No stall |   1.185|
| YML058W   | SML1      | No stall |   1.177|
| YDR339C   | FCF1      | No stall |   1.174|
| YGL254W   | FZF1      | No stall |   1.172|
| YPL252C   | YAH1      | No stall |   1.172|
| YNL124W   | NAF1      | Stall    |   1.166|
| YDL223C   | HBT1      | No stall |   1.165|
| YML022W   | APT1      | No stall |   1.155|
| YKL222C   | NA        | Stall    |   1.149|
| YHR122W   | CIA2      | No stall |   1.139|
| YCL044C   | MGR1      | No stall |   1.138|
| YDL012C   | NA        | No stall |   1.138|
| YEL046C   | GLY1      | No stall |   1.137|
| YLR059C   | REX2      | No stall |   1.134|
| YJL158C   | CIS3      | No stall |   1.131|
| YBL082C   | ALG3      | No stall |   1.123|
| YAR050W   | FLO1      | No stall |   1.119|
| YNR036C   | MRPS12    | No stall |   1.117|
| YCL036W   | GFD2      | No stall |   1.115|
| YEL048C   | TCA17     | No stall |   1.114|
| YER088C   | DOT6      | No stall |   1.112|
| YHR049W   | FSH1      | No stall |   1.107|
| YML042W   | CAT2      | No stall |   1.107|
| YOL154W   | ZPS1      | No stall |   1.107|
| YDL014W   | NOP1      | No stall |   1.098|
| YKL198C   | PTK1      | Stall    |   1.096|
| YBL016W   | FUS3      | No stall |   1.086|
| YLR375W   | STP3      | Stall    |   1.083|
| YBR031W   | RPL4A     | No stall |   1.069|
| YPL163C   | SVS1      | No stall |   1.069|
| YLR058C   | SHM2      | No stall |   1.065|
| YNL141W   | AAH1      | No stall |   1.060|
| YLR439W   | MRPL4     | No stall |   1.059|
| YBL054W   | TOD6      | No stall |   1.052|
| YCL004W   | PGS1      | No stall |   1.048|
| YLL051C   | FRE6      | No stall |   1.045|
| YGR280C   | PXR1      | Stall    |   1.044|
| YHR148W   | IMP3      | No stall |   1.044|
| YML132W   | COS3      | No stall |   1.040|
| YJR004C   | SAG1      | No stall |   1.035|
| YEL019C   | MMS21     | No stall |   1.033|
| YMR016C   | SOK2      | No stall |   1.031|
| YKL084W   | HOT13     | No stall |   1.025|
| YGR156W   | PTI1      | Stall    |   1.022|
| YBR231C   | SWC5      | No stall |   1.011|
| YBL028C   | NA        | No stall |   1.011|
| YLR388W   | RPS29A    | No stall |   1.010|
| YEL047C   | FRD1      | No stall |   1.006|
| YBR173C   | UMP1      | No stall |   1.006|
| YOL119C   | MCH4      | No stall |   1.005|
| YHR211W   | FLO5      | No stall |   1.000|
| YCR072C   | RSA4      | No stall |   0.998|
| YER083C   | GET2      | Stall    |   0.997|
| YER109C   | FLO8      | No stall |   0.994|
| YGL028C   | SCW11     | No stall |   0.992|
| YGL222C   | EDC1      | No stall |   0.992|
| YPR144C   | NOC4      | No stall |   0.987|
| YHL040C   | ARN1      | No stall |   0.982|
| YMR180C   | CTL1      | No stall |   0.981|
| YGL172W   | NUP49     | No stall |   0.978|
| YDR006C   | SOK1      | No stall |   0.977|
| YMR058W   | FET3      | No stall |   0.972|
| YJR161C   | COS5      | No stall |   0.971|
| YMR149W   | SWP1      | No stall |   0.966|
| YLL019C   | KNS1      | Stall    |   0.966|
| YNL071W   | LAT1      | No stall |   0.962|
| YOL158C   | ENB1      | No stall |   0.960|
| YDR101C   | ARX1      | No stall |   0.958|
| YGR159C   | NSR1      | No stall |   0.958|
| YOR329C   | SCD5      | Stall    |   0.955|
| YGL156W   | AMS1      | No stall |   0.952|
| YOR120W   | GCY1      | No stall |   0.949|
| YPL089C   | RLM1      | No stall |   0.948|
| YLR237W   | THI7      | No stall |   0.944|
| YCR065W   | HCM1      | No stall |   0.944|
| YIR023W   | DAL81     | No stall |   0.943|
| YCR107W   | AAD3      | No stall |   0.943|
| YOR303W   | CPA1      | No stall |   0.942|
| YNL336W   | COS1      | No stall |   0.942|
| YIL117C   | PRM5      | No stall |   0.939|
| YMR318C   | ADH6      | No stall |   0.939|
| YNL025C   | SSN8      | No stall |   0.937|
| YDL168W   | SFA1      | No stall |   0.934|
| YGL178W   | MPT5      | No stall |   0.933|
| YKL078W   | DHR2      | No stall |   0.930|
| YDR059C   | UBC5      | No stall |   0.930|
| YKL004W   | AUR1      | No stall |   0.925|
| YGL056C   | SDS23     | No stall |   0.925|
| YLR403W   | SFP1      | No stall |   0.925|
| YHL007C   | STE20     | Stall    |   0.924|
| YPL036W   | PMA2      | No stall |   0.914|
| YHR094C   | HXT1      | No stall |   0.910|
| YIL123W   | SIM1      | No stall |   0.909|
| YFL034C-B | MOB2      | Stall    |   0.909|
| YGR191W   | HIP1      | No stall |   0.908|
| YDR012W   | RPL4B     | No stall |   0.907|
| YOR359W   | VTS1      | No stall |   0.901|
| YOL028C   | YAP7      | Stall    |   0.900|
| YIL046W   | MET30     | No stall |   0.896|
| YGR212W   | SLI1      | No stall |   0.896|
| YDR398W   | UTP5      | Stall    |   0.893|
| YEL009C   | GCN4      | Stall    |   0.886|
| YJL041W   | NSP1      | No stall |   0.886|
| YKL032C   | IXR1      | Stall    |   0.883|
| YDL189W   | RBS1      | Stall    |   0.881|
| YOL067C   | RTG1      | Stall    |   0.881|
| YBL075C   | SSA3      | No stall |   0.880|
| YDR420W   | HKR1      | Stall    |   0.878|
| YLR092W   | SUL2      | No stall |   0.876|
| YHR135C   | YCK1      | No stall |   0.875|
| YCR048W   | ARE1      | No stall |   0.872|
| YOR251C   | TUM1      | No stall |   0.867|
| YOL104C   | NDJ1      | No stall |   0.867|
| YPL128C   | TBF1      | No stall |   0.863|
| YKR050W   | TRK2      | Stall    |   0.860|
| YFR036W   | CDC26     | No stall |   0.858|
| YEL043W   | NA        | No stall |   0.857|
| YJR119C   | JHD2      | No stall |   0.852|
| YOR349W   | CIN1      | No stall |   0.852|
| YHR208W   | BAT1      | No stall |   0.851|
| YJL025W   | RRN7      | No stall |   0.849|
| YNL145W   | MFA2      | No stall |   0.842|
| YER152C   | NA        | No stall |   0.842|
| YHR206W   | SKN7      | No stall |   0.841|
| YOL128C   | YGK3      | No stall |   0.840|
| YNL310C   | ZIM17     | No stall |   0.840|
| YDR207C   | UME6      | Stall    |   0.836|
| YHR199C   | AIM46     | No stall |   0.832|
| YLR370C   | ARC18     | No stall |   0.832|
| YHR021C   | RPS27B    | No stall |   0.831|
| YER009W   | NTF2      | No stall |   0.831|
| YHR179W   | OYE2      | No stall |   0.831|
| YPR085C   | ASA1      | No stall |   0.830|
| YLR074C   | BUD20     | Stall    |   0.826|
| YDL052C   | SLC1      | No stall |   0.826|
| YDL220C   | CDC13     | No stall |   0.824|
| YPL204W   | HRR25     | No stall |   0.823|
| YMR070W   | MOT3      | No stall |   0.821|
| YLR437C   | DIF1      | No stall |   0.818|
| YJL084C   | ALY2      | No stall |   0.815|
| YGR195W   | SKI6      | No stall |   0.814|
| YGL166W   | CUP2      | Stall    |   0.813|
| YOR315W   | SFG1      | No stall |   0.812|
| YDR260C   | SWM1      | No stall |   0.812|
| YOR140W   | SFL1      | No stall |   0.812|
| YDL088C   | ASM4      | No stall |   0.809|
| YKR046C   | PET10     | No stall |   0.809|
| YJL151C   | SNA3      | No stall |   0.808|
| YCR059C   | YIH1      | No stall |   0.808|
| YKR055W   | RHO4      | No stall |   0.807|
| YLR401C   | DUS3      | Stall    |   0.805|
| YGR234W   | YHB1      | No stall |   0.800|
| YLR431C   | ATG23     | No stall |   0.799|
| YGR173W   | RBG2      | No stall |   0.798|
| YNL283C   | WSC2      | No stall |   0.795|
| YGL071W   | AFT1      | Stall    |   0.794|
| YDL087C   | LUC7      | Stall    |   0.792|
| YPR159W   | KRE6      | No stall |   0.792|
| YMR099C   | NA        | No stall |   0.792|
| YML027W   | YOX1      | Stall    |   0.792|
| YIL003W   | CFD1      | No stall |   0.791|
| YJR063W   | RPA12     | No stall |   0.790|
| YPL177C   | CUP9      | No stall |   0.788|
| YNL158W   | PGA1      | No stall |   0.788|
| YLR047C   | FRE8      | No stall |   0.787|
| YCL008C   | STP22     | Stall    |   0.785|
| YCL029C   | BIK1      | No stall |   0.785|
| YOR023C   | AHC1      | Stall    |   0.781|
| YCL058W-A | ADF1      | No stall |   0.781|
| YML113W   | DAT1      | Stall    |   0.776|
| YLR110C   | CCW12     | No stall |   0.775|
| YKL019W   | RAM2      | No stall |   0.774|
| Q0130     | OLI1      | No stall |   0.773|
| YKL082C   | RRP14     | Stall    |   0.772|
| YAL040C   | CLN3      | No stall |   0.771|
| YLL013C   | PUF3      | Stall    |   0.771|
| YGL234W   | ADE5,7    | No stall |   0.771|
| YHR205W   | SCH9      | No stall |   0.770|
| YLR027C   | AAT2      | No stall |   0.769|
| YPR094W   | RDS3      | No stall |   0.769|
| YML114C   | TAF8      | No stall |   0.769|
| YBR034C   | HMT1      | No stall |   0.768|
| YJL065C   | DLS1      | No stall |   0.766|
| YLR394W   | CST9      | Stall    |   0.766|
| YOL043C   | NTG2      | Stall    |   0.766|
| YCR067C   | SED4      | Stall    |   0.764|
| YMR266W   | RSN1      | Stall    |   0.763|
| YDL243C   | AAD4      | No stall |   0.762|
| YMR164C   | MSS11     | Stall    |   0.762|
| YLR328W   | NMA1      | Stall    |   0.761|
| YJL121C   | RPE1      | No stall |   0.761|
| YIL133C   | RPL16A    | No stall |   0.757|
| YML014W   | TRM9      | No stall |   0.754|
| YLR157C   | ASP3-1    | No stall |   0.751|
| YER116C   | SLX8      | Stall    |   0.751|
| YOR381W   | FRE3      | No stall |   0.750|
| YDR267C   | CIA1      | No stall |   0.749|
| YAR015W   | ADE1      | No stall |   0.747|
| YPL169C   | MEX67     | No stall |   0.746|
| YGL211W   | NCS6      | No stall |   0.745|
| YGR146C   | ECL1      | No stall |   0.742|
| YOR355W   | GDS1      | No stall |   0.739|
| YCL009C   | ILV6      | No stall |   0.737|
| YKL164C   | PIR1      | No stall |   0.736|
| YHR071W   | PCL5      | No stall |   0.734|
| YNR012W   | URK1      | No stall |   0.733|
| YNL230C   | ELA1      | Stall    |   0.732|
| YDR363W   | ESC2      | No stall |   0.727|
| YDL197C   | ASF2      | No stall |   0.727|
| YMR037C   | MSN2      | Stall    |   0.727|
| YOR278W   | HEM4      | No stall |   0.727|
| YCL039W   | GID7      | No stall |   0.727|
| YDR232W   | HEM1      | No stall |   0.725|
| YLR008C   | PAM18     | Stall    |   0.724|
| YNL075W   | IMP4      | No stall |   0.722|
| YDR513W   | GRX2      | No stall |   0.722|
| YGL251C   | HFM1      | Stall    |   0.718|
| YOL007C   | CSI2      | No stall |   0.715|
| YLR131C   | ACE2      | No stall |   0.715|
| YKL051W   | SFK1      | No stall |   0.715|
| YAL012W   | CYS3      | No stall |   0.713|
| YDR451C   | YHP1      | Stall    |   0.713|
| YPL016W   | SWI1      | Stall    |   0.712|
| YBR238C   | NA        | No stall |   0.711|
| YKL068W   | NUP100    | No stall |   0.711|
| YNL014W   | HEF3      | Stall    |   0.706|
| YOR028C   | CIN5      | No stall |   0.705|
| YBL066C   | SEF1      | No stall |   0.705|
| YHR040W   | BCD1      | No stall |   0.703|
| YPL076W   | GPI2      | No stall |   0.702|
| YGL215W   | CLG1      | No stall |   0.701|
| YOR307C   | SLY41     | No stall |   0.700|
| YGR014W   | MSB2      | No stall |   0.700|
| YER180C   | ISC10     | No stall |   0.698|
| YHR127W   | NA        | No stall |   0.697|
| YLR120C   | YPS1      | No stall |   0.696|
| YOR148C   | SPP2      | Stall    |   0.696|
| YER154W   | OXA1      | No stall |   0.696|
| YHR088W   | RPF1      | Stall    |   0.696|
| YLR155C   | ASP3-1    | No stall |   0.695|
| YML107C   | PML39     | No stall |   0.695|
| YDL080C   | THI3      | No stall |   0.694|
| YJR041C   | URB2      | No stall |   0.693|
| YIL135C   | VHS2      | No stall |   0.693|
| YBL014C   | RRN6      | Stall    |   0.690|
| YHR072W   | ERG7      | No stall |   0.688|
| YBR293W   | VBA2      | No stall |   0.688|
| YDR171W   | HSP42     | No stall |   0.687|
| YML051W   | GAL80     | No stall |   0.687|
| YER132C   | PMD1      | Stall    |   0.687|
| YIL145C   | PAN6      | No stall |   0.685|
| YLR009W   | RLP24     | Stall    |   0.684|
| YLR390W-A | CCW14     | No stall |   0.684|
| YOR188W   | MSB1      | No stall |   0.682|
| YER062C   | GPP2      | No stall |   0.682|
| YER086W   | ILV1      | No stall |   0.680|
| YHR199C-A | NBL1      | Stall    |   0.679|
| YOL155C   | HPF1      | No stall |   0.679|
| YAL036C   | RBG1      | No stall |   0.678|
| YDR441C   | APT2      | No stall |   0.676|
| YOR316C   | COT1      | No stall |   0.676|
| YFL024C   | EPL1      | Stall    |   0.675|
| YER040W   | GLN3      | Stall    |   0.673|
| YGL083W   | SCY1      | No stall |   0.673|
| YDL224C   | WHI4      | No stall |   0.672|
| YKR094C   | RPL40B    | No stall |   0.672|
| YEL030W   | ECM10     | No stall |   0.672|
| YGL014W   | PUF4      | Stall    |   0.672|
| YBL051C   | PIN4      | Stall    |   0.668|
| YOR208W   | PTP2      | No stall |   0.667|
| YDR077W   | SED1      | No stall |   0.665|
| YHR102W   | KIC1      | No stall |   0.665|
| YBR302C   | COS3      | No stall |   0.664|
| YNL053W   | MSG5      | No stall |   0.664|
| YLR116W   | MSL5      | No stall |   0.664|
| YDR533C   | HSP31     | No stall |   0.663|
| YML005W   | TRM12     | No stall |   0.662|
| YLR006C   | SSK1      | No stall |   0.661|
| YML128C   | MSC1      | No stall |   0.660|
| YGR158C   | MTR3      | No stall |   0.660|
| YDR411C   | DFM1      | No stall |   0.659|
| YPR181C   | SEC23     | No stall |   0.659|
| YHL033C   | RPL8A     | No stall |   0.658|
| YLR055C   | SPT8      | No stall |   0.656|
| YHR061C   | GIC1      | No stall |   0.656|
| YKL192C   | ACP1      | No stall |   0.655|
| YEL066W   | HPA3      | No stall |   0.654|
| YGR032W   | GSC2      | Stall    |   0.654|
| YGR060W   | ERG25     | No stall |   0.653|
| YER052C   | HOM3      | No stall |   0.653|
| YPR093C   | ASR1      | No stall |   0.653|
| YPR198W   | SGE1      | No stall |   0.653|
| YGR197C   | SNG1      | No stall |   0.652|
| YDR461W   | MFA1      | No stall |   0.652|
| YOL141W   | PPM2      | Stall    |   0.651|
| YGL122C   | NAB2      | No stall |   0.651|
| YNL197C   | WHI3      | No stall |   0.651|
| YPL258C   | THI21     | No stall |   0.651|
| YLR318W   | EST2      | No stall |   0.650|
| YAL053W   | FLC2      | No stall |   0.649|
| YMR065W   | KAR5      | No stall |   0.647|
| YOR299W   | BUD7      | No stall |   0.646|
| YBR167C   | POP7      | No stall |   0.646|
| YGL032C   | AGA2      | No stall |   0.644|
| YBR261C   | TAE1      | No stall |   0.644|
| YAL020C   | ATS1      | No stall |   0.644|
| YFR022W   | ROG3      | No stall |   0.641|
| YBR035C   | PDX3      | No stall |   0.640|
| YNL062C   | GCD10     | No stall |   0.639|
| YNL066W   | SUN4      | No stall |   0.639|
| YLR199C   | PBA1      | No stall |   0.637|
| YLR332W   | MID2      | No stall |   0.637|
| YOR204W   | DED1      | No stall |   0.636|
| YMR002W   | MIX17     | No stall |   0.636|
| YCL028W   | RNQ1      | No stall |   0.634|
| YJL083W   | TAX4      | Stall    |   0.633|
| YOL072W   | THP1      | No stall |   0.632|
| YKL055C   | OAR1      | No stall |   0.631|
| YJL138C   | TIF2      | Stall    |   0.630|
| YKR071C   | DRE2      | Stall    |   0.630|
| YEL044W   | IES6      | No stall |   0.627|
| YNR038W   | DBP6      | No stall |   0.627|
| YLR034C   | SMF3      | No stall |   0.627|
| YGR049W   | SCM4      | No stall |   0.626|
| YLR139C   | SLS1      | No stall |   0.624|
| YBR061C   | TRM7      | No stall |   0.621|
| YBR212W   | NGR1      | No stall |   0.621|
| YJR091C   | JSN1      | No stall |   0.618|
| YER065C   | ICL1      | No stall |   0.617|
| YHR074W   | QNS1      | No stall |   0.616|
| YIR031C   | DAL7      | No stall |   0.616|
| YBR215W   | HPC2      | No stall |   0.616|
| YIL087C   | AIM19     | No stall |   0.616|
| YNL069C   | RPL16B    | No stall |   0.614|
| YLR160C   | ASP3-1    | No stall |   0.613|
| YLR176C   | RFX1      | No stall |   0.612|
| YFL017W-A | SMX2      | No stall |   0.611|
| YOR210W   | RPB10     | No stall |   0.611|
| YLR052W   | IES3      | Stall    |   0.611|
| YML059C   | NTE1      | No stall |   0.611|
| YBL063W   | KIP1      | No stall |   0.609|
| YPL061W   | ALD6      | No stall |   0.609|
| YNL044W   | YIP3      | No stall |   0.609|
| YOL116W   | MSN1      | No stall |   0.608|
| YOR272W   | YTM1      | Stall    |   0.607|
| YIL148W   | RPL40B    | No stall |   0.607|
| YKL054C   | DEF1      | No stall |   0.606|
| YGL169W   | SUA5      | No stall |   0.606|
| YOL130W   | ALR1      | No stall |   0.605|
| YOR066W   | MSA1      | No stall |   0.604|
| YJR155W   | AAD10     | Stall    |   0.603|
| YGL126W   | SCS3      | No stall |   0.603|
| YDR036C   | EHD3      | No stall |   0.602|
| YLR096W   | KIN2      | No stall |   0.601|
| YEL053C   | MAK10     | No stall |   0.600|
| YDR472W   | TRS31     | No stall |   0.600|
| YOR163W   | DDP1      | No stall |   0.598|
| YMR079W   | SEC14     | No stall |   0.598|
| YLR256W   | HAP1      | Stall    |   0.598|
| YLR091W   | GEP5      | No stall |   0.596|
| YOL091W   | SPO21     | No stall |   0.595|
| YCR028C   | FEN2      | No stall |   0.594|
| YLL035W   | GRC3      | No stall |   0.594|
| YML101C   | CUE4      | No stall |   0.592|
| YPL179W   | PPQ1      | Stall    |   0.592|
| YNR046W   | TRM112    | No stall |   0.591|
| YHR039C   | MSC7      | No stall |   0.591|
| YBL074C   | AAR2      | No stall |   0.591|
| YJR021C   | REC107    | Stall    |   0.591|
| YOL061W   | PRS5      | No stall |   0.590|
| YBR133C   | HSL7      | No stall |   0.590|
| YNL076W   | MKS1      | No stall |   0.589|
| YLR121C   | YPS3      | No stall |   0.587|
| YJL109C   | UTP10     | Stall    |   0.586|
| YER127W   | LCP5      | Stall    |   0.586|
| YPL164C   | MLH3      | No stall |   0.586|
| YDR384C   | ATO3      | No stall |   0.585|
| YER002W   | NOP16     | Stall    |   0.585|
| YFL044C   | OTU1      | No stall |   0.584|
| YNL309W   | STB1      | No stall |   0.581|
| YPL126W   | NAN1      | Stall    |   0.581|
| YML112W   | CTK3      | No stall |   0.580|
| YDL203C   | ACK1      | No stall |   0.579|
| YNL047C   | SLM2      | No stall |   0.578|
| YLR457C   | NBP1      | Stall    |   0.578|
| YMR068W   | AVO2      | No stall |   0.578|
| YKL098W   | MTC2      | No stall |   0.578|
| YOR091W   | TMA46     | Stall    |   0.578|
| YNL245C   | CWC25     | No stall |   0.576|
| YOR081C   | TGL5      | No stall |   0.575|
| YER169W   | RPH1      | Stall    |   0.575|
| YLR150W   | STM1      | No stall |   0.575|
| YER010C   | NA        | No stall |   0.574|
| YKR061W   | KTR2      | No stall |   0.573|
| YER041W   | YEN1      | Stall    |   0.573|
| YJR097W   | JJJ3      | No stall |   0.571|
| YKL038W   | RGT1      | Stall    |   0.570|
| YHR084W   | STE12     | No stall |   0.570|
| YNL327W   | EGT2      | No stall |   0.569|
| YGR010W   | NMA2      | No stall |   0.569|
| YER004W   | FMP52     | No stall |   0.569|
| YIL019W   | FAF1      | Stall    |   0.568|
| YGL162W   | SUT1      | No stall |   0.567|
| YDR103W   | STE5      | No stall |   0.567|
| YER038C   | KRE29     | No stall |   0.566|
| YMR255W   | GFD1      | No stall |   0.565|
| YNL186W   | UBP10     | Stall    |   0.565|
| YHR042W   | NCP1      | No stall |   0.565|
| YKL029C   | MAE1      | No stall |   0.564|
| YML052W   | SUR7      | No stall |   0.564|
| YGR216C   | GPI1      | No stall |   0.562|
| YIL140W   | AXL2      | No stall |   0.561|
| YNL314W   | DAL82     | No stall |   0.559|
| YOL164W   | BDS1      | No stall |   0.559|
| YOR337W   | TEA1      | No stall |   0.559|
| YKL043W   | PHD1      | No stall |   0.557|
| YDR409W   | SIZ1      | No stall |   0.556|
| YDR496C   | PUF6      | No stall |   0.555|
| YJL162C   | JJJ2      | Stall    |   0.555|
| YLR355C   | ILV5      | No stall |   0.554|
| YLR397C   | AFG2      | No stall |   0.553|
| YJL058C   | BIT61     | No stall |   0.551|
| YHR032W   | ERC1      | No stall |   0.550|
| YBR158W   | AMN1      | No stall |   0.550|
| YMR306W   | FKS3      | No stall |   0.550|
| YJL069C   | UTP18     | No stall |   0.549|
| YER144C   | UBP5      | No stall |   0.547|
| YDR309C   | GIC2      | No stall |   0.547|
| YEL040W   | UTR2      | No stall |   0.546|
| YNL251C   | NRD1      | Stall    |   0.546|
| YOR328W   | PDR10     | Stall    |   0.546|
| YDR021W   | FAL1      | No stall |   0.545|
| YML095C   | RAD10     | No stall |   0.544|
| YML019W   | OST6      | No stall |   0.544|
| YIL047C   | SYG1      | Stall    |   0.543|
| YGR081C   | SLX9      | Stall    |   0.543|
| YPL031C   | PHO85     | No stall |   0.542|
| YGL231C   | EMC4      | No stall |   0.542|
| YKR060W   | UTP30     | No stall |   0.542|
| YGR097W   | ASK10     | No stall |   0.541|
| YNL182C   | IPI3      | No stall |   0.540|
| YHL031C   | GOS1      | No stall |   0.540|
| YMR038C   | CCS1      | No stall |   0.540|
| YER115C   | SPR6      | Stall    |   0.538|
| YNR048W   | NA        | No stall |   0.538|
| YOL068C   | HST1      | Stall    |   0.537|
| YLR158C   | ASP3-1    | No stall |   0.536|
| YMR170C   | ALD2      | No stall |   0.536|
| YNL112W   | DBP2      | No stall |   0.535|
| YNL110C   | NOP15     | No stall |   0.535|
| YDL060W   | TSR1      | No stall |   0.534|
| YJL062W   | LAS21     | No stall |   0.533|
| YAR042W   | SWH1      | No stall |   0.533|
| YCR057C   | PWP2      | No stall |   0.532|
| YDR489W   | SLD5      | No stall |   0.531|
| YOR063W   | RPL3      | No stall |   0.530|
| YKL159C   | RCN1      | No stall |   0.529|
| YER060W   | FCY21     | No stall |   0.528|
| YMR064W   | AEP1      | No stall |   0.527|
| YNR013C   | PHO91     | No stall |   0.526|
| YDR165W   | TRM82     | No stall |   0.525|
| YJL100W   | LSB6      | No stall |   0.525|
| YPL037C   | EGD1      | No stall |   0.525|
| YGR057C   | LST7      | No stall |   0.525|
| YOR128C   | ADE2      | No stall |   0.523|
| YMR303C   | ADH2      | No stall |   0.522|
| YOL008W   | COQ10     | No stall |   0.522|
| YIL156W   | UBP7      | Stall    |   0.522|
| YKL012W   | PRP40     | No stall |   0.522|
| YLR231C   | BNA5      | No stall |   0.521|
| YMR300C   | ADE4      | No stall |   0.521|
| YMR042W   | ARG80     | No stall |   0.519|
| YML081C-A | ATP18     | No stall |   0.519|
| YLR412W   | BER1      | Stall    |   0.519|
| YPL127C   | HHO1      | No stall |   0.519|
| YMR093W   | UTP15     | Stall    |   0.518|
| YER114C   | BOI2      | Stall    |   0.518|
| YDR082W   | STN1      | No stall |   0.517|
| YEL003W   | GIM4      | No stall |   0.517|
| YKR025W   | RPC37     | No stall |   0.517|
| YJL033W   | HCA4      | Stall    |   0.517|
| YMR297W   | PRC1      | No stall |   0.516|
| YDR233C   | RTN1      | No stall |   0.516|
| YHR169W   | DBP8      | Stall    |   0.515|
| YMR217W   | GUA1      | No stall |   0.513|
| YLR192C   | HCR1      | Stall    |   0.512|
| YHR193C   | EGD2      | No stall |   0.512|
| YMR290C   | HAS1      | No stall |   0.512|
| YGL018C   | JAC1      | No stall |   0.512|
| YDL179W   | PCL9      | No stall |   0.511|
| YKR004C   | ECM9      | No stall |   0.510|
| YPR041W   | TIF5      | No stall |   0.510|
| YJR016C   | ILV3      | No stall |   0.510|
| YOL093W   | TRM10     | Stall    |   0.509|
| YDR075W   | PPH3      | No stall |   0.509|
| YML099C   | ARG81     | Stall    |   0.509|
| YDR032C   | PST2      | No stall |   0.509|
| YPL184C   | MRN1      | No stall |   0.508|
| YBR267W   | REI1      | Stall    |   0.508|
| YLR215C   | CDC123    | No stall |   0.507|
| YBR295W   | PCA1      | No stall |   0.507|
| YPR171W   | BSP1      | Stall    |   0.507|
| YMR276W   | DSK2      | No stall |   0.507|
| YCR079W   | PTC6      | No stall |   0.507|
| YIL120W   | QDR1      | No stall |   0.507|
| YGL073W   | HSF1      | No stall |   0.506|
| YDR064W   | RPS13     | No stall |   0.505|
| YLR228C   | ECM22     | No stall |   0.505|
| YPL032C   | SVL3      | Stall    |   0.504|
| YNL138W-A | YSF3      | No stall |   0.504|
| YIL136W   | OM45      | No stall |   0.504|
| YPL129W   | TAF14     | No stall |   0.502|
| YBR214W   | SDS24     | No stall |   0.502|
| YCR060W   | TAH1      | No stall |   0.502|
| YKR076W   | ECM4      | No stall |   0.502|
| YBR156C   | SLI15     | Stall    |   0.501|
| YJL118W   | NA        | Stall    |   0.501|
| YOL100W   | PKH2      | Stall    |   0.500|
| YBR193C   | MED8      | No stall |   0.500|
| YGL090W   | LIF1      | Stall    |   0.500|
| YNL094W   | APP1      | Stall    |   0.500|
| YCR038C   | BUD5      | No stall |   0.500|
| YGR211W   | ZPR1      | No stall |   0.499|
| YJL133W   | MRS3      | No stall |   0.498|
| YGL030W   | RPL30     | No stall |   0.497|
| YLR206W   | ENT2      | Stall    |   0.497|
| YDR285W   | ZIP1      | Stall    |   0.497|
| YPR065W   | ROX1      | Stall    |   0.497|
| YNL063W   | MTQ1      | No stall |   0.497|
| YDR011W   | SNQ2      | Stall    |   0.495|
| YBR030W   | RKM3      | No stall |   0.495|
| YNL147W   | LSM7      | No stall |   0.495|
| YOR224C   | RPB8      | No stall |   0.494|
| YPR016C   | TIF6      | No stall |   0.493|
| YJR148W   | BAT2      | No stall |   0.492|
| YGL113W   | SLD3      | Stall    |   0.492|
| YDR404C   | RPB7      | No stall |   0.492|
| YLL027W   | ISA1      | Stall    |   0.492|
| YGR251W   | NOP19     | Stall    |   0.491|
| YBL092W   | RPL32     | Stall    |   0.491|
| YJR135C   | MCM22     | No stall |   0.491|
| YGL256W   | ADH4      | No stall |   0.490|
| YOL124C   | TRM11     | No stall |   0.490|
| YLR106C   | MDN1      | Stall    |   0.490|
| YBR112C   | CYC8      | No stall |   0.489|
| YGR172C   | YIP1      | No stall |   0.488|
| YLL008W   | DRS1      | Stall    |   0.488|
| YOR002W   | ALG6      | No stall |   0.488|
| YNL006W   | LST8      | No stall |   0.488|
| YBR223C   | TDP1      | No stall |   0.487|
| YBL085W   | BOI1      | No stall |   0.487|
| YGR128C   | UTP8      | No stall |   0.487|
| YGL174W   | BUD13     | Stall    |   0.486|
| YCR091W   | KIN82     | Stall    |   0.486|
| YCR084C   | TUP1      | No stall |   0.486|
| YOL051W   | GAL11     | No stall |   0.485|
| YJR120W   | NA        | No stall |   0.485|
| YCR020C   | PET18     | No stall |   0.484|
| YGR133W   | PEX4      | No stall |   0.484|
| YOR011W   | AUS1      | No stall |   0.484|
| YAL043C   | PTA1      | No stall |   0.484|
| YLR205C   | HMX1      | No stall |   0.483|
| YLR316C   | TAD3      | No stall |   0.483|
| YHR161C   | YAP1801   | No stall |   0.483|
| YDR242W   | AMD2      | No stall |   0.483|
| YNL078W   | NIS1      | Stall    |   0.482|
| YIR001C   | SGN1      | No stall |   0.481|
| YBR264C   | YPT10     | No stall |   0.481|
| YPR184W   | GDB1      | No stall |   0.480|
| YPL193W   | RSA1      | Stall    |   0.480|
| YIL130W   | ASG1      | No stall |   0.480|
| YNL221C   | POP1      | Stall    |   0.478|
| YFL026W   | STE2      | No stall |   0.477|
| YDL248W   | COS7      | No stall |   0.476|
| YJR106W   | ECM27     | No stall |   0.475|
| YER151C   | UBP3      | Stall    |   0.475|
| YGR006W   | PRP18     | No stall |   0.474|
| YBL039C   | URA7      | No stall |   0.474|
| YOR229W   | WTM2      | No stall |   0.473|
| YIR034C   | LYS1      | No stall |   0.473|
| YEL062W   | NPR2      | No stall |   0.472|
| YOR197W   | MCA1      | Stall    |   0.472|
| YHR072W-A | NOP10     | No stall |   0.471|
| YPL008W   | CHL1      | Stall    |   0.471|
| YAL021C   | CCR4      | No stall |   0.471|
| YNL123W   | NMA111    | No stall |   0.471|
| YNL097C   | PHO23     | Stall    |   0.470|
| YPL085W   | SEC16     | Stall    |   0.470|
| YDR456W   | NHX1      | No stall |   0.469|
| YNL241C   | ZWF1      | No stall |   0.468|
| YBR082C   | UBC4      | No stall |   0.468|
| YNL152W   | INN1      | No stall |   0.468|
| YDR155C   | CPR1      | No stall |   0.466|
| YMR008C   | PLB1      | No stall |   0.466|
| YML063W   | RPS1B     | Stall    |   0.466|
| YER035W   | EDC2      | Stall    |   0.465|
| YOL113W   | SKM1      | No stall |   0.465|
| YOL030W   | GAS5      | Stall    |   0.465|
| YFL022C   | FRS2      | No stall |   0.465|
| YBR108W   | AIM3      | Stall    |   0.465|
| YER145C   | FTR1      | No stall |   0.464|
| YOR231W   | MKK1      | No stall |   0.464|
| YCR088W   | ABP1      | Stall    |   0.463|
| YDR359C   | EAF1      | No stall |   0.463|
| YOR221C   | MCT1      | No stall |   0.462|
| YGL057C   | GEP7      | No stall |   0.461|
| YML043C   | RRN11     | Stall    |   0.461|
| YBL060W   | YEL1      | No stall |   0.460|
| YIR026C   | YVH1      | No stall |   0.460|
| YOR252W   | TMA16     | Stall    |   0.460|
| YOR004W   | UTP23     | Stall    |   0.459|
| YPL112C   | PEX25     | No stall |   0.459|
| YJR046W   | TAH11     | No stall |   0.459|
| YDL167C   | NRP1      | No stall |   0.458|
| YHR020W   | NA        | No stall |   0.458|
| YGR185C   | TYS1      | Stall    |   0.458|
| YHL011C   | PRS3      | No stall |   0.456|
| YNL298W   | CLA4      | No stall |   0.456|
| YGL078C   | DBP3      | Stall    |   0.455|
| YDR186C   | SND1      | No stall |   0.452|
| YOR145C   | PNO1      | No stall |   0.451|
| YDL128W   | VCX1      | No stall |   0.450|
| YHR089C   | GAR1      | Stall    |   0.450|
| YDL033C   | SLM3      | No stall |   0.450|
| YJR072C   | NPA3      | Stall    |   0.450|
| YAL038W   | CDC19     | No stall |   0.450|
| YBR029C   | CDS1      | No stall |   0.449|
| YPR070W   | MED1      | No stall |   0.449|
| YOR227W   | HER1      | Stall    |   0.449|
| YNL103W   | MET4      | Stall    |   0.448|
| YIL030C   | SSM4      | No stall |   0.447|
| YOR056C   | NOB1      | Stall    |   0.447|
| YHL047C   | ARN2      | No stall |   0.446|
| YKL014C   | URB1      | No stall |   0.446|
| YAL014C   | SYN8      | No stall |   0.446|
| YJL191W   | RPS14B    | Stall    |   0.445|
| YDR310C   | SUM1      | Stall    |   0.445|
| YNL288W   | CAF40     | No stall |   0.445|
| YOR267C   | HRK1      | Stall    |   0.444|
| YGR056W   | RSC1      | Stall    |   0.444|
| YDR263C   | DIN7      | Stall    |   0.443|
| YOR030W   | DFG16     | Stall    |   0.443|
| YLR044C   | PDC1      | No stall |   0.443|
| YHR030C   | SLT2      | No stall |   0.443|
| YKL127W   | PGM1      | No stall |   0.442|
| YLR451W   | LEU3      | Stall    |   0.442|
| YLR026C   | SED5      | No stall |   0.442|
| YMR033W   | ARP9      | No stall |   0.441|
| YDR120C   | TRM1      | Stall    |   0.441|
| YFL021W   | GAT1      | Stall    |   0.441|
| YHR111W   | UBA4      | No stall |   0.439|
| YOL031C   | SIL1      | No stall |   0.439|
| YKL020C   | SPT23     | Stall    |   0.438|
| YCR106W   | RDS1      | Stall    |   0.438|
| YGL044C   | RNA15     | No stall |   0.437|
| YDR206W   | EBS1      | No stall |   0.437|
| YNL055C   | POR1      | No stall |   0.434|
| YNL052W   | COX5A     | No stall |   0.434|
| YPR112C   | MRD1      | Stall    |   0.434|
| YGL167C   | PMR1      | Stall    |   0.434|
| YCL055W   | KAR4      | No stall |   0.433|
| YGR090W   | UTP22     | No stall |   0.433|
| YLR209C   | PNP1      | No stall |   0.433|
| YML017W   | PSP2      | No stall |   0.433|
| YKL090W   | CUE2      | No stall |   0.433|
| YDR140W   | MTQ2      | No stall |   0.433|
| YOL052C   | SPE2      | No stall |   0.432|
| YLR060W   | FRS1      | No stall |   0.432|
| YDR160W   | SSY1      | No stall |   0.432|
| YDR390C   | UBA2      | No stall |   0.432|
| YKR021W   | ALY1      | Stall    |   0.431|
| YER082C   | UTP7      | No stall |   0.431|
| YKL146W   | AVT3      | No stall |   0.431|
| YOR319W   | HSH49     | No stall |   0.431|
| YJR093C   | FIP1      | Stall    |   0.430|
| YNL175C   | NOP13     | Stall    |   0.429|
| YGR077C   | PEX8      | No stall |   0.429|
| YGR080W   | TWF1      | Stall    |   0.429|
| YDR205W   | MSC2      | No stall |   0.428|
| YJL149W   | DAS1      | No stall |   0.428|
| YMR049C   | ERB1      | No stall |   0.427|
| YGR283C   | NA        | Stall    |   0.427|
| YDR272W   | GLO2      | No stall |   0.427|
| YNL334C   | SNO2      | No stall |   0.426|
| YIL020C   | HIS6      | No stall |   0.425|
| YLR172C   | DPH5      | No stall |   0.424|
| YNL119W   | NCS2      | No stall |   0.424|
| YBL025W   | RRN10     | No stall |   0.423|
| YLR459W   | GAB1      | No stall |   0.421|
| YPR115W   | RGC1      | No stall |   0.421|
| YNL219C   | ALG9      | No stall |   0.420|
| YDL150W   | RPC53     | Stall    |   0.420|
| YGL120C   | PRP43     | Stall    |   0.420|
| YMR305C   | SCW10     | No stall |   0.420|
| YBR067C   | TIP1      | No stall |   0.420|
| YEL036C   | ANP1      | No stall |   0.419|
| YFR046C   | CNN1      | No stall |   0.419|
| YDR200C   | VPS64     | No stall |   0.419|
| YOR161C   | PNS1      | Stall    |   0.419|
| YPL246C   | RBD2      | No stall |   0.419|
| YNL282W   | POP3      | No stall |   0.419|
| YJL090C   | DPB11     | No stall |   0.419|
| YOR373W   | NUD1      | No stall |   0.418|
| YDR087C   | RRP1      | No stall |   0.417|
| YKR092C   | SRP40     | No stall |   0.417|
| YOR098C   | NUP1      | No stall |   0.417|
| YDR172W   | SUP35     | No stall |   0.416|
| YJL134W   | LCB3      | No stall |   0.416|
| YMR242C   | RPL20A    | No stall |   0.415|
| YNR015W   | SMM1      | No stall |   0.414|
| YOR071C   | NRT1      | No stall |   0.414|
| YEL063C   | CAN1      | No stall |   0.414|
| YDL213C   | NOP6      | Stall    |   0.412|
| YPR042C   | PUF2      | Stall    |   0.412|
| YOR322C   | LDB19     | No stall |   0.412|
| YBL064C   | PRX1      | No stall |   0.411|
| YER110C   | KAP123    | No stall |   0.411|
| YLR373C   | VID22     | No stall |   0.410|
| YER087C-B | SBH1      | No stall |   0.410|
| YOL055C   | THI20     | No stall |   0.409|
| YIR006C   | PAN1      | Stall    |   0.409|
| YHR006W   | STP2      | Stall    |   0.409|
| YJL117W   | PHO86     | Stall    |   0.409|
| YPR040W   | TIP41     | No stall |   0.408|
| YNL031C   | HHT2      | No stall |   0.408|
| YLR056W   | ERG3      | No stall |   0.407|
| YOR380W   | RDR1      | Stall    |   0.407|
| YDR144C   | MKC7      | No stall |   0.407|
| YLR214W   | FRE1      | No stall |   0.407|
| YNL238W   | KEX2      | Stall    |   0.406|
| YOL125W   | TRM13     | Stall    |   0.406|
| YGR142W   | BTN2      | No stall |   0.404|
| YER068W   | MOT2      | Stall    |   0.404|
| YDL048C   | STP4      | Stall    |   0.404|
| YLR147C   | SMD3      | No stall |   0.404|
| YGL139W   | FLC3      | No stall |   0.403|
| YER019W   | ISC1      | No stall |   0.403|
| YNR075W   | COS10     | No stall |   0.403|
| YGR058W   | PEF1      | No stall |   0.402|
| YDR287W   | INM2      | No stall |   0.402|
| YCR069W   | CPR4      | No stall |   0.402|
| YLL018C   | DPS1      | No stall |   0.401|
| YNR059W   | MNT4      | No stall |   0.401|
| YHR197W   | RIX1      | No stall |   0.400|
| YPR106W   | ISR1      | No stall |   0.400|
| YGL025C   | PGD1      | Stall    |   0.399|
| YCR004C   | YCP4      | No stall |   0.399|
| YNR054C   | ESF2      | No stall |   0.398|
| YOR202W   | HIS3      | No stall |   0.397|
| YJL179W   | PFD1      | No stall |   0.397|
| YBR055C   | PRP6      | Stall    |   0.397|
| YBR191W   | RPL21A    | No stall |   0.397|
| YJL174W   | KRE9      | No stall |   0.394|
| YLR017W   | MEU1      | No stall |   0.394|
| YKR044W   | UIP5      | Stall    |   0.394|
| YLL012W   | YEH1      | No stall |   0.393|
| YDR153C   | ENT5      | Stall    |   0.392|
| YDR091C   | RLI1      | No stall |   0.392|
| YLR084C   | RAX2      | Stall    |   0.391|
| YDR152W   | GIR2      | No stall |   0.391|
| YMR113W   | FOL3      | No stall |   0.390|
| YOR003W   | YSP3      | No stall |   0.390|
| YPL093W   | NOG1      | No stall |   0.390|
| YOR046C   | DBP5      | No stall |   0.389|
| YDL138W   | RGT2      | Stall    |   0.389|
| YHR116W   | COX23     | No stall |   0.389|
| YOL020W   | TAT2      | No stall |   0.389|
| YPR010C   | RPA135    | No stall |   0.389|
| YPL202C   | AFT2      | No stall |   0.389|
| YOR226C   | ISU2      | No stall |   0.389|
| YPL217C   | BMS1      | Stall    |   0.389|
| YMR296C   | LCB1      | No stall |   0.389|
| YER153C   | PET122    | No stall |   0.388|
| YAL048C   | GEM1      | Stall    |   0.388|
| YGL243W   | TAD1      | No stall |   0.387|
| YKL021C   | MAK11     | Stall    |   0.387|
| YHR064C   | SSZ1      | No stall |   0.386|
| YGL136C   | MRM2      | No stall |   0.386|
| YER126C   | NSA2      | Stall    |   0.385|
| YIL154C   | IMP2'     | Stall    |   0.385|
| YLR094C   | GIS3      | Stall    |   0.385|
| YJL001W   | PRE3      | No stall |   0.385|
| YIL023C   | YKE4      | No stall |   0.385|
| YHR204W   | MNL1      | No stall |   0.384|
| YMR106C   | YKU80     | No stall |   0.383|
| YOR253W   | NAT5      | No stall |   0.383|
| YLR452C   | SST2      | No stall |   0.383|
| YPL038W   | MET31     | No stall |   0.382|
| YNL249C   | MPA43     | No stall |   0.382|
| YAL022C   | FUN26     | No stall |   0.382|
| YIL031W   | ULP2      | No stall |   0.382|
| YDL103C   | QRI1      | No stall |   0.381|
| YPR137W   | RRP9      | Stall    |   0.381|
| YLR196W   | PWP1      | Stall    |   0.381|
| YLR186W   | EMG1      | No stall |   0.381|
| YER100W   | UBC6      | No stall |   0.381|
| YDR483W   | KRE2      | No stall |   0.380|
| YDL201W   | TRM8      | Stall    |   0.380|
| YLR449W   | FPR4      | Stall    |   0.380|
| YNL256W   | FOL1      | No stall |   0.380|
| YKR026C   | GCN3      | No stall |   0.379|
| YHR083W   | SAM35     | No stall |   0.378|
| YOR360C   | PDE2      | No stall |   0.378|
| YAL039C   | CYC3      | No stall |   0.377|
| YHR149C   | SKG6      | No stall |   0.377|
| YNL308C   | KRI1      | Stall    |   0.377|
| YOR043W   | WHI2      | Stall    |   0.377|
| YOR377W   | ATF1      | No stall |   0.377|
| YOL138C   | RTC1      | No stall |   0.377|
| YOR181W   | LAS17     | Stall    |   0.376|
| YOR162C   | YRR1      | No stall |   0.376|
| YIR012W   | SQT1      | No stall |   0.376|
| YPR164W   | MMS1      | No stall |   0.375|
| YOR020C   | HSP10     | No stall |   0.375|
| YNR043W   | MVD1      | No stall |   0.375|
| YNL207W   | RIO2      | Stall    |   0.374|
| YJR086W   | STE18     | No stall |   0.374|
| YPL091W   | GLR1      | No stall |   0.374|
| YLR043C   | TRX1      | No stall |   0.374|
| YBL033C   | RIB1      | No stall |   0.374|
| YAL047C   | SPC72     | No stall |   0.374|
| YOR239W   | ABP140    | Stall    |   0.373|
| YDR361C   | BCP1      | Stall    |   0.372|
| YDR379W   | RGA2      | No stall |   0.372|
| YOR233W   | KIN4      | Stall    |   0.372|
| YDR463W   | STP1      | Stall    |   0.371|
| YHL035C   | VMR1      | No stall |   0.370|
| YBR036C   | CSG2      | No stall |   0.370|
| YPL159C   | PET20     | Stall    |   0.369|
| YBR249C   | ARO4      | No stall |   0.368|
| YGL004C   | RPN14     | No stall |   0.368|
| YKR043C   | SHB17     | No stall |   0.367|
| YOR007C   | SGT2      | No stall |   0.367|
| YMR319C   | FET4      | No stall |   0.367|
| YNL074C   | MLF3      | No stall |   0.367|
| YOL060C   | MAM3      | No stall |   0.367|
| YAL025C   | MAK16     | Stall    |   0.366|
| YGL253W   | HXK2      | No stall |   0.366|
| YPR017C   | DSS4      | No stall |   0.365|
| YJR109C   | CPA2      | No stall |   0.365|
| YNL331C   | AAD14     | No stall |   0.364|
| YDR110W   | FOB1      | Stall    |   0.364|
| YMR067C   | UBX4      | No stall |   0.364|
| YGR170W   | PSD2      | Stall    |   0.362|
| YHR013C   | ARD1      | No stall |   0.362|
| YOL081W   | IRA2      | No stall |   0.362|
| YPL166W   | ATG29     | No stall |   0.362|
| YBR265W   | TSC10     | No stall |   0.362|
| YER129W   | SAK1      | No stall |   0.361|
| YHR100C   | GEP4      | No stall |   0.361|
| YCR081W   | SRB8      | No stall |   0.361|
| YGL070C   | RPB9      | No stall |   0.360|
| YER081W   | SER3      | No stall |   0.360|
| YGR208W   | SER2      | No stall |   0.360|
| YBR286W   | APE3      | No stall |   0.360|
| YGL009C   | LEU1      | No stall |   0.360|
| YKL017C   | HCS1      | Stall    |   0.360|
| YCR073W-A | SOL2      | Stall    |   0.359|
| YKR077W   | MSA2      | No stall |   0.359|
| YPR030W   | CSR2      | Stall    |   0.358|
| YGL043W   | DST1      | No stall |   0.357|
| YGL098W   | USE1      | No stall |   0.356|
| YEL042W   | GDA1      | No stall |   0.356|
| YGR264C   | MES1      | No stall |   0.355|
| YDR224C   | HTB1      | Stall    |   0.355|
| YCL063W   | VAC17     | No stall |   0.355|
| YAL010C   | MDM10     | No stall |   0.354|
| YNR006W   | VPS27     | Stall    |   0.354|
| YOR008C   | SLG1      | No stall |   0.354|
| YDR297W   | SUR2      | No stall |   0.354|
| YAL019W   | FUN30     | No stall |   0.354|
| YBR252W   | DUT1      | No stall |   0.353|
| YDR422C   | SIP1      | No stall |   0.353|
| YFR050C   | PRE4      | No stall |   0.353|
| YGR037C   | ACB1      | No stall |   0.353|
| YIL061C   | SNP1      | No stall |   0.353|
| YDR098C   | GRX3      | No stall |   0.351|
| YJR132W   | NMD5      | No stall |   0.351|
| YKL184W   | SPE1      | No stall |   0.351|
| YLR227C   | ADY4      | Stall    |   0.350|
| YBR106W   | PHO88     | No stall |   0.349|
| YNL248C   | RPA49     | No stall |   0.349|
| YAR019C   | CDC15     | No stall |   0.349|
| YPL055C   | LGE1      | No stall |   0.348|
| YNL054W   | VAC7      | No stall |   0.347|
| YGL062W   | PYC1      | Stall    |   0.346|
| YOR363C   | PIP2      | No stall |   0.345|
| YKL180W   | RPL17A    | No stall |   0.344|
| YJR001W   | AVT1      | No stall |   0.344|
| YGR177C   | ATF2      | No stall |   0.343|
| YNL042W   | BOP3      | Stall    |   0.343|
| YGR083C   | GCD2      | Stall    |   0.343|
| YGL103W   | RPL28     | No stall |   0.343|
| YHR010W   | RPL27A    | Stall    |   0.342|
| YKL174C   | TPO5      | No stall |   0.342|
| YGL154C   | LYS5      | No stall |   0.342|
| YJR007W   | SUI2      | No stall |   0.342|
| YIL164C   | NIT1      | No stall |   0.341|
| YNL167C   | SKO1      | Stall    |   0.341|
| YMR312W   | ELP6      | No stall |   0.340|
| YKL081W   | TEF4      | No stall |   0.340|
| YFL017C   | GNA1      | No stall |   0.340|
| YDR214W   | AHA1      | No stall |   0.340|
| YLR090W   | XDJ1      | No stall |   0.339|
| YBR086C   | IST2      | No stall |   0.339|
| YLR144C   | ACF2      | No stall |   0.338|
| YER176W   | ECM32     | Stall    |   0.338|
| YJR082C   | EAF6      | No stall |   0.338|
| YDR143C   | SAN1      | No stall |   0.338|
| YER093C   | TSC11     | No stall |   0.338|
| YKR024C   | DBP7      | Stall    |   0.338|
| YDR492W   | IZH1      | No stall |   0.338|
| YPL012W   | RRP12     | Stall    |   0.337|
| YGL225W   | VRG4      | No stall |   0.337|
| YLR208W   | SEC13     | No stall |   0.336|
| YLR083C   | EMP70     | No stall |   0.336|
| YMR080C   | NAM7      | No stall |   0.336|
| YGR227W   | DIE2      | Stall    |   0.335|
| YNR018W   | RCF2      | No stall |   0.335|
| YHR200W   | RPN10     | No stall |   0.335|
| YDL230W   | PTP1      | No stall |   0.334|
| YBR155W   | CNS1      | No stall |   0.334|
| YDR228C   | PCF11     | No stall |   0.334|
| YCL034W   | LSB5      | No stall |   0.334|
| YMR015C   | ERG5      | No stall |   0.334|
| YBR101C   | FES1      | No stall |   0.334|
| YLR133W   | CKI1      | No stall |   0.334|
| YBR079C   | RPG1      | Stall    |   0.334|
| YDL111C   | RRP42     | No stall |   0.333|
| YKL176C   | LST4      | No stall |   0.332|
| YOR108W   | LEU9      | No stall |   0.332|
| YDL235C   | YPD1      | No stall |   0.332|
| YDR123C   | INO2      | No stall |   0.331|
| YCR077C   | PAT1      | Stall    |   0.331|
| YMR172W   | HOT1      | Stall    |   0.331|
| YER112W   | LSM4      | No stall |   0.331|
| YPR121W   | THI22     | No stall |   0.331|
| YOL133W   | HRT1      | Stall    |   0.331|
| YLR418C   | CDC73     | No stall |   0.330|
| YGL197W   | MDS3      | Stall    |   0.329|
| YKL188C   | PXA2      | No stall |   0.329|
| YDR423C   | CAD1      | Stall    |   0.329|
| YMR236W   | TAF9      | No stall |   0.329|
| YMR153W   | NUP53     | No stall |   0.328|
| YGL190C   | CDC55     | No stall |   0.327|
| YER056C   | FCY2      | No stall |   0.327|
| YOL144W   | NOP8      | Stall    |   0.327|
| YOR288C   | MPD1      | No stall |   0.327|
| YDR443C   | SSN2      | No stall |   0.327|
| YPL131W   | RPL5      | No stall |   0.326|
| YEL005C   | VAB2      | No stall |   0.326|
| YJR067C   | YAE1      | No stall |   0.324|
| YGR074W   | SMD1      | No stall |   0.324|
| YJL056C   | ZAP1      | No stall |   0.324|
| YFR021W   | ATG18     | No stall |   0.324|
| YBR152W   | SPP381    | No stall |   0.324|
| YPR032W   | SRO7      | No stall |   0.323|
| YHR046C   | INM1      | No stall |   0.323|
| YDR365C   | ESF1      | Stall    |   0.323|
| YER107C   | GLE2      | No stall |   0.323|
| YLR114C   | AVL9      | No stall |   0.322|
| YBR289W   | SNF5      | Stall    |   0.322|
| YMR123W   | PKR1      | No stall |   0.322|
| YOR006C   | TSR3      | No stall |   0.322|
| YDL030W   | PRP9      | Stall    |   0.322|
| YMR019W   | STB4      | Stall    |   0.322|
| YPR034W   | ARP7      | No stall |   0.321|
| YDR299W   | BFR2      | No stall |   0.320|
| YJL023C   | PET130    | No stall |   0.319|
| YDR092W   | UBC13     | No stall |   0.319|
| YJR083C   | ACF4      | No stall |   0.319|
| YAL032C   | PRP45     | No stall |   0.319|
| YDR539W   | FDC1      | No stall |   0.318|
| YMR215W   | GAS3      | No stall |   0.318|
| YBR263W   | SHM1      | No stall |   0.317|
| YJR002W   | MPP10     | Stall    |   0.317|
| YOR358W   | HAP5      | No stall |   0.317|
| YKL103C   | APE1      | No stall |   0.317|
| YHL021C   | AIM17     | No stall |   0.317|
| YLR134W   | PDC5      | No stall |   0.316|
| YBL015W   | ACH1      | No stall |   0.316|
| YFL033C   | RIM15     | Stall    |   0.316|
| YDR436W   | PPZ2      | No stall |   0.316|
| YKL008C   | LAC1      | Stall    |   0.316|
| YNL068C   | FKH2      | No stall |   0.316|
| YDL005C   | MED2      | No stall |   0.315|
| YDR328C   | SKP1      | No stall |   0.315|
| YGR162W   | TIF4631   | Stall    |   0.314|
| YBL011W   | SCT1      | Stall    |   0.314|
| YDR465C   | RMT2      | No stall |   0.314|
| YPL157W   | TGS1      | Stall    |   0.314|
| YGL099W   | LSG1      | Stall    |   0.313|
| YDL208W   | NHP2      | No stall |   0.313|
| YLR359W   | ADE13     | No stall |   0.313|
| YHR115C   | DMA1      | No stall |   0.312|
| YLR347C   | KAP95     | No stall |   0.312|
| YLR350W   | ORM2      | No stall |   0.312|
| YOR340C   | RPA43     | No stall |   0.311|
| YDR395W   | SXM1      | No stall |   0.311|
| YPL137C   | GIP3      | No stall |   0.311|
| YMR043W   | MCM1      | No stall |   0.311|
| YLR014C   | PPR1      | Stall    |   0.311|
| YBR188C   | NTC20     | No stall |   0.311|
| YPL043W   | NOP4      | Stall    |   0.310|
| YPR009W   | SUT2      | Stall    |   0.310|
| YAR027W   | UIP3      | No stall |   0.310|
| YNL153C   | GIM3      | No stall |   0.310|
| YDR121W   | DPB4      | No stall |   0.310|
| YGL236C   | MTO1      | No stall |   0.310|
| YMR222C   | FSH2      | No stall |   0.309|
| YIL049W   | DFG10     | No stall |   0.309|
| YGL053W   | PRM8      | No stall |   0.307|
| YER167W   | BCK2      | Stall    |   0.307|
| YKR065C   | PAM17     | No stall |   0.307|
| YER178W   | PDA1      | No stall |   0.306|
| YCR039C   | HMLALPHA2 | No stall |   0.306|
| YMR298W   | LIP1      | No stall |   0.306|
| YLR384C   | IKI3      | Stall    |   0.306|
| YDL217C   | TIM22     | No stall |   0.305|
| YNL104C   | LEU4      | No stall |   0.305|
| YBR198C   | TAF5      | No stall |   0.304|
| YCL010C   | SGF29     | No stall |   0.304|
| YMR243C   | ZRC1      | No stall |   0.304|
| YLR018C   | POM34     | No stall |   0.303|
| YNL261W   | ORC5      | Stall    |   0.303|
| YIL142W   | CCT2      | No stall |   0.302|
| YGR113W   | DAM1      | No stall |   0.302|
| YPL183C   | RTT10     | No stall |   0.302|
| YHR056C   | NA        | No stall |   0.302|
| YPR079W   | MRL1      | No stall |   0.302|
| YGR200C   | ELP2      | No stall |   0.302|
| YDR034C   | LYS14     | Stall    |   0.301|
| YOL010W   | RCL1      | No stall |   0.301|
| YDL051W   | LHP1      | Stall    |   0.301|
| YPR168W   | NUT2      | No stall |   0.301|
| YLR376C   | PSY3      | No stall |   0.301|
| YHR163W   | SOL3      | No stall |   0.301|
| YJL122W   | ALB1      | Stall    |   0.301|
| YOL111C   | MDY2      | No stall |   0.300|
| YDR406W   | PDR15     | No stall |   0.300|
| YBR253W   | SRB6      | No stall |   0.300|
| YOR123C   | LEO1      | Stall    |   0.300|
| YOR269W   | PAC1      | No stall |   0.299|
| YGR260W   | TNA1      | No stall |   0.299|
| YPR190C   | RPC82     | No stall |   0.299|
| YKL099C   | UTP11     | Stall    |   0.299|
| YGR245C   | SDA1      | Stall    |   0.299|
| YLR327C   | TMA10     | No stall |   0.299|
| YGL144C   | ROG1      | No stall |   0.298|
| YOR301W   | RAX1      | No stall |   0.298|
| YMR146C   | TIF34     | No stall |   0.298|
| YDR270W   | CCC2      | No stall |   0.298|
| YNL113W   | RPC19     | No stall |   0.298|
| YGR038W   | ORM1      | No stall |   0.298|
| YLR436C   | ECM30     | No stall |   0.297|
| YHL003C   | LAG1      | No stall |   0.297|
| YBR291C   | CTP1      | No stall |   0.297|
| YPL180W   | TCO89     | No stall |   0.297|
| YNR026C   | SEC12     | No stall |   0.297|
| YDR293C   | SSD1      | Stall    |   0.296|
| YOR027W   | STI1      | No stall |   0.296|
| YML075C   | HMG1      | No stall |   0.296|
| YJL048C   | UBX6      | No stall |   0.296|
| YFR044C   | DUG1      | No stall |   0.296|
| YKL211C   | TRP3      | No stall |   0.296|
| YGR181W   | TIM13     | No stall |   0.295|
| YJL011C   | RPC17     | No stall |   0.294|
| YIR008C   | PRI1      | No stall |   0.294|
| YPL117C   | IDI1      | No stall |   0.293|
| YHR052W   | CIC1      | Stall    |   0.293|
| YDR324C   | UTP4      | Stall    |   0.293|
| YGL087C   | MMS2      | Stall    |   0.292|
| YKR099W   | BAS1      | Stall    |   0.292|
| YGR106C   | VOA1      | No stall |   0.292|
| YJL153C   | INO1      | No stall |   0.291|
| YMR234W   | RNH1      | No stall |   0.290|
| YKR081C   | RPF2      | No stall |   0.290|
| YBR109C   | CMD1      | No stall |   0.290|
| YER122C   | GLO3      | No stall |   0.289|
| YMR277W   | FCP1      | No stall |   0.289|
| YDR189W   | SLY1      | No stall |   0.289|
| YBR107C   | IML3      | No stall |   0.289|
| YOR312C   | RPL20A    | No stall |   0.289|
| YDR221W   | GTB1      | Stall    |   0.289|
| YLR425W   | TUS1      | Stall    |   0.289|
| YDR538W   | PAD1      | No stall |   0.288|
| YCR034W   | ELO2      | No stall |   0.288|
| YNL236W   | SIN4      | No stall |   0.288|
| YAL028W   | FRT2      | No stall |   0.287|
| YLR219W   | MSC3      | Stall    |   0.287|
| YFL050C   | ALR2      | No stall |   0.287|
| YLR381W   | CTF3      | No stall |   0.287|
| YJR096W   | NA        | No stall |   0.286|
| YIL033C   | BCY1      | No stall |   0.286|
| YOR167C   | RPS28A    | No stall |   0.286|
| YPR122W   | AXL1      | Stall    |   0.286|
| YDR202C   | RAV2      | No stall |   0.285|
| YOR361C   | PRT1      | Stall    |   0.285|
| YPR074C   | TKL1      | No stall |   0.284|
| YJL014W   | CCT3      | No stall |   0.284|
| YFL002C   | SPB4      | Stall    |   0.284|
| YGR233C   | PHO81     | No stall |   0.284|
| YKL183W   | LOT5      | No stall |   0.284|
| YNL275W   | BOR1      | No stall |   0.283|
| YLR293C   | GSP1      | No stall |   0.282|
| YBR194W   | AIM4      | No stall |   0.281|
| YLR229C   | CDC42     | No stall |   0.281|
| YKL191W   | DPH2      | No stall |   0.281|
| YOR305W   | RRG7      | No stall |   0.281|
| YGL210W   | YPT32     | No stall |   0.281|
| YOR189W   | IES4      | No stall |   0.280|
| YKR042W   | UTH1      | No stall |   0.280|
| YNR049C   | MSO1      | No stall |   0.280|
| YJL172W   | CPS1      | No stall |   0.280|
| YMR203W   | TOM40     | No stall |   0.280|
| YML006C   | GIS4      | No stall |   0.279|
| YNL111C   | CYB5      | No stall |   0.279|
| YDR188W   | CCT6      | Stall    |   0.279|
| YGR214W   | RPS0A     | No stall |   0.277|
| YOL112W   | MSB4      | Stall    |   0.277|
| YDR499W   | LCD1      | No stall |   0.277|
| YMR273C   | ZDS1      | Stall    |   0.276|
| YJL201W   | ECM25     | No stall |   0.276|
| YKR027W   | BCH2      | Stall    |   0.276|
| YNL180C   | RHO5      | Stall    |   0.275|
| YER124C   | DSE1      | Stall    |   0.275|
| YML130C   | ERO1      | No stall |   0.275|
| YMR039C   | SUB1      | Stall    |   0.274|
| YDR195W   | REF2      | No stall |   0.274|
| YLR443W   | ECM7      | No stall |   0.274|
| YIR033W   | MGA2      | Stall    |   0.274|
| YNL164C   | IBD2      | No stall |   0.273|
| YMR239C   | RNT1      | No stall |   0.273|
| YIL147C   | SLN1      | No stall |   0.273|
| YJR105W   | ADO1      | No stall |   0.273|
| YDR078C   | SHU2      | No stall |   0.272|
| YER123W   | YCK3      | No stall |   0.271|
| YOR243C   | PUS7      | No stall |   0.271|
| YPL190C   | NAB3      | Stall    |   0.271|
| YMR131C   | RRB1      | No stall |   0.271|
| YNL132W   | KRE33     | Stall    |   0.271|
| YDR128W   | MTC5      | No stall |   0.271|
| YMR167W   | MLH1      | No stall |   0.271|
| YBR172C   | SMY2      | Stall    |   0.270|
| YHL023C   | NPR3      | Stall    |   0.270|
| YEL032W   | MCM3      | Stall    |   0.270|
| YOR112W   | CEX1      | No stall |   0.270|
| YPR026W   | ATH1      | No stall |   0.270|
| YKL106W   | AAT1      | No stall |   0.269|
| YGR124W   | ASN2      | No stall |   0.269|
| YKL072W   | STB6      | Stall    |   0.269|
| YJL212C   | OPT1      | No stall |   0.268|
| YLR129W   | DIP2      | Stall    |   0.268|
| YIL091C   | UTP25     | Stall    |   0.267|
| YCR047C   | BUD23     | Stall    |   0.267|
| YDR480W   | DIG2      | No stall |   0.267|
| YDR530C   | APA2      | No stall |   0.267|
| YPL160W   | CDC60     | No stall |   0.266|
| YMR201C   | RAD14     | Stall    |   0.266|
| YHR155W   | LAM1      | No stall |   0.265|
| YOR058C   | ASE1      | No stall |   0.265|
| YPL214C   | THI6      | No stall |   0.263|
| YER090W   | TRP2      | No stall |   0.263|
| YOR109W   | INP53     | Stall    |   0.263|
| YKL156W   | RPS27A    | No stall |   0.263|
| YLR020C   | YEH2      | Stall    |   0.263|
| YDL099W   | BUG1      | Stall    |   0.262|
| YIL149C   | MLP2      | No stall |   0.262|
| YPL213W   | LEA1      | No stall |   0.262|
| YMR110C   | HFD1      | No stall |   0.261|
| YNL329C   | PEX6      | No stall |   0.261|
| YHR187W   | IKI1      | No stall |   0.261|
| YCR008W   | SAT4      | No stall |   0.260|
| YBR210W   | ERV15     | No stall |   0.260|
| YPL065W   | VPS28     | No stall |   0.260|
| YNL023C   | FAP1      | No stall |   0.260|
| YML102W   | CAC2      | No stall |   0.259|
| YDL013W   | SLX5      | Stall    |   0.259|
| YLL023C   | POM33     | No stall |   0.259|
| YNL149C   | PGA2      | Stall    |   0.259|
| YLR146C   | SPE4      | No stall |   0.258|
| YJL050W   | MTR4      | Stall    |   0.258|
| YHR069C   | RRP4      | No stall |   0.258|
| YLR361C   | DCR2      | No stall |   0.258|
| YLR002C   | NOC3      | Stall    |   0.258|
| YBR078W   | ECM33     | No stall |   0.258|
| YNR027W   | BUD17     | No stall |   0.258|
| YGL189C   | RPS26A    | No stall |   0.257|
| YPR020W   | ATP20     | No stall |   0.257|
| YIL109C   | SEC24     | Stall    |   0.257|
| YOL077C   | BRX1      | No stall |   0.256|
| YOR219C   | STE13     | No stall |   0.256|
| YKR072C   | SIS2      | No stall |   0.256|
| YLR432W   | IMD3      | No stall |   0.256|
| YPR161C   | SGV1      | No stall |   0.256|
| YER033C   | ZRG8      | Stall    |   0.255|
| YDR208W   | MSS4      | Stall    |   0.255|
| YDR268W   | MSW1      | No stall |   0.255|
| YOR179C   | SYC1      | No stall |   0.254|
| YKR091W   | SRL3      | No stall |   0.254|
| YHR029C   | YHI9      | No stall |   0.254|
| YLL045C   | RPL8B     | No stall |   0.254|
| YBR004C   | GPI18     | No stall |   0.254|
| YER047C   | SAP1      | Stall    |   0.254|
| YPR072W   | NOT5      | Stall    |   0.253|
| YBL005W   | PDR3      | Stall    |   0.253|
| YAL058W   | CNE1      | No stall |   0.253|
| YDR040C   | ENA1      | No stall |   0.253|
| YBR163W   | EXO5      | No stall |   0.253|
| YOR040W   | GLO4      | No stall |   0.253|
| YHR066W   | SSF1      | Stall    |   0.252|
| YNL199C   | GCR2      | No stall |   0.252|
| YNL016W   | PUB1      | No stall |   0.252|
| YDL124W   | NA        | No stall |   0.252|
| YKL005C   | BYE1      | Stall    |   0.251|
| YDL017W   | CDC7      | No stall |   0.251|
| YBR222C   | PCS60     | No stall |   0.250|
| YLR265C   | NEJ1      | Stall    |   0.250|
| YIL052C   | RPL34B    | No stall |   0.250|
| YPL253C   | VIK1      | No stall |   0.250|
| YDL153C   | SAS10     | Stall    |   0.250|
| YNL116W   | DMA2      | No stall |   0.249|
| YML093W   | UTP14     | Stall    |   0.249|
| YMR308C   | PSE1      | No stall |   0.249|
| YER054C   | GIP2      | No stall |   0.249|
| YJR062C   | NTA1      | No stall |   0.249|
| YBR247C   | ENP1      | No stall |   0.249|
| YMR285C   | NGL2      | Stall    |   0.248|
| YGR145W   | ENP2      | Stall    |   0.248|
| YGR155W   | CYS4      | No stall |   0.248|
| YBR276C   | PPS1      | No stall |   0.247|
| YGL023C   | PIB2      | No stall |   0.247|
| YJL010C   | NOP9      | Stall    |   0.247|
| YNL292W   | PUS4      | Stall    |   0.247|
| YLL026W   | HSP104    | No stall |   0.247|
| YIL103W   | DPH1      | No stall |   0.247|
| YDR386W   | MUS81     | Stall    |   0.247|
| YDR351W   | SBE2      | No stall |   0.247|
| YDL195W   | SEC31     | No stall |   0.246|
| YKL052C   | ASK1      | No stall |   0.246|
| YHR086W   | NAM8      | No stall |   0.246|
| YER165W   | PAB1      | No stall |   0.246|
| YHR101C   | BIG1      | No stall |   0.246|
| YKL128C   | PMU1      | No stall |   0.246|
| YNL036W   | NCE103    | No stall |   0.246|
| YDR425W   | SNX41     | No stall |   0.246|
| YGR078C   | PAC10     | No stall |   0.246|
| YIL062C   | ARC15     | No stall |   0.246|
| YLL009C   | COX17     | No stall |   0.246|
| YOR206W   | NOC2      | Stall    |   0.245|
| YNR010W   | CSE2      | No stall |   0.245|
| YDR192C   | NUP42     | No stall |   0.245|
| YFR005C   | SAD1      | No stall |   0.245|
| YHR175W   | CTR2      | No stall |   0.244|
| YML035C   | AMD1      | No stall |   0.244|
| YKL135C   | APL2      | No stall |   0.244|
| YER057C   | HMF1      | No stall |   0.244|
| YEL037C   | RAD23     | No stall |   0.244|
| YOL059W   | GPD2      | No stall |   0.244|
| YHR008C   | SOD2      | No stall |   0.243|
| YDR028C   | REG1      | Stall    |   0.243|
| YJL019W   | MPS3      | No stall |   0.242|
| YJL080C   | SCP160    | No stall |   0.242|
| YNL154C   | YCK2      | No stall |   0.241|
| YOR341W   | RPA190    | Stall    |   0.241|
| YBR288C   | APM3      | No stall |   0.241|
| YFL045C   | SEC53     | No stall |   0.241|
| YPL273W   | SAM4      | No stall |   0.241|
| YML103C   | NUP188    | No stall |   0.241|
| YKR059W   | TIF2      | Stall    |   0.241|
| YBR026C   | ETR1      | No stall |   0.241|
| YGR198W   | YPP1      | No stall |   0.241|
| YNL317W   | PFS2      | No stall |   0.240|
| YKL025C   | PAN3      | Stall    |   0.240|
| YLR135W   | SLX4      | Stall    |   0.240|
| YBR142W   | MAK5      | Stall    |   0.239|
| YNL339C   | YRF1-6    | Stall    |   0.239|
| YBL056W   | PTC3      | No stall |   0.239|
| YKL094W   | YJU3      | No stall |   0.239|
| YIL143C   | SSL2      | Stall    |   0.238|
| YNL079C   | TPM1      | No stall |   0.238|
| YPL018W   | CTF19     | No stall |   0.238|
| YLL038C   | ENT4      | No stall |   0.238|
| YJL141C   | YAK1      | No stall |   0.237|
| YFR051C   | RET2      | Stall    |   0.237|
| YIL110W   | HPM1      | No stall |   0.237|
| YKR080W   | MTD1      | No stall |   0.237|
| YGL181W   | GTS1      | No stall |   0.237|
| YKL095W   | YJU2      | No stall |   0.236|
| YGL237C   | HAP2      | No stall |   0.236|
| YDR174W   | HMO1      | Stall    |   0.235|
| YNL163C   | RIA1      | No stall |   0.235|
| YPL234C   | VMA11     | No stall |   0.234|
| YBL026W   | LSM2      | No stall |   0.234|
| YAL002W   | VPS8      | No stall |   0.233|
| YBR119W   | MUD1      | Stall    |   0.232|
| YDR227W   | SIR4      | No stall |   0.232|
| YPR143W   | RRP15     | Stall    |   0.232|
| YML092C   | PRE8      | No stall |   0.232|
| YDR354W   | TRP4      | No stall |   0.231|
| YPL105C   | SYH1      | Stall    |   0.230|
| YGR068C   | ART5      | No stall |   0.230|
| YBR195C   | MSI1      | No stall |   0.230|
| YDR367W   | KEI1      | No stall |   0.230|
| YBL071W-A | KTI11     | No stall |   0.230|
| YDR410C   | STE14     | No stall |   0.229|
| YER075C   | PTP3      | No stall |   0.229|
| YKR038C   | KAE1      | No stall |   0.228|
| YPR119W   | CLB2      | No stall |   0.228|
| YPL086C   | ELP3      | No stall |   0.227|
| YGL049C   | TIF4632   | Stall    |   0.227|
| YNR003C   | RPC34     | No stall |   0.227|
| YJL210W   | PEX2      | No stall |   0.226|
| YDL231C   | BRE4      | No stall |   0.226|
| YJL127C   | SPT10     | Stall    |   0.225|
| YLR342W   | FKS1      | Stall    |   0.224|
| YOR205C   | GEP3      | No stall |   0.224|
| YHR062C   | RPP1      | No stall |   0.224|
| YBL031W   | SHE1      | Stall    |   0.224|
| YOL135C   | MED7      | No stall |   0.224|
| YPL049C   | DIG1      | No stall |   0.224|
| YNL254C   | RTC4      | No stall |   0.224|
| YDL063C   | SYO1      | No stall |   0.223|
| YGL171W   | ROK1      | Stall    |   0.223|
| YER161C   | SPT2      | Stall    |   0.223|
| YGL244W   | RTF1      | Stall    |   0.223|
| YNR035C   | ARC35     | No stall |   0.223|
| YDR412W   | RRP17     | Stall    |   0.223|
| YAR002W   | NUP60     | No stall |   0.222|
| YDR247W   | VHS1      | No stall |   0.222|
| YOR038C   | HIR2      | No stall |   0.222|
| YDR501W   | PLM2      | Stall    |   0.222|
| YOL140W   | ARG8      | No stall |   0.221|
| YGL039W   | NA        | No stall |   0.221|
| YER159C   | BUR6      | No stall |   0.221|
| YBR105C   | VID24     | No stall |   0.221|
| YBL084C   | CDC27     | No stall |   0.220|
| YLR354C   | TAL1      | No stall |   0.220|
| YNL007C   | SIS1      | Stall    |   0.220|
| YOL009C   | MDM12     | No stall |   0.220|
| YHL002W   | HSE1      | No stall |   0.219|
| YMR135C   | GID8      | No stall |   0.219|
| YHR099W   | TRA1      | Stall    |   0.219|
| YDL215C   | GDH2      | No stall |   0.219|
| YGL163C   | RAD54     | No stall |   0.218|
| YOR294W   | RRS1      | Stall    |   0.217|
| YOR103C   | OST2      | No stall |   0.217|
| YDR326C   | YSP2      | Stall    |   0.216|
| YBL093C   | ROX3      | No stall |   0.215|
| YLL050C   | COF1      | No stall |   0.215|
| YOR293W   | RPS10A    | No stall |   0.215|
| YHR058C   | MED6      | No stall |   0.215|
| YHL025W   | SNF6      | No stall |   0.215|
| YLR233C   | EST1      | No stall |   0.215|
| YEL031W   | SPF1      | No stall |   0.213|
| YHR134W   | WSS1      | No stall |   0.212|
| YNL263C   | YIF1      | No stall |   0.212|
| YOR119C   | RIO1      | Stall    |   0.212|
| YDR226W   | ADK1      | No stall |   0.212|
| YPL040C   | ISM1      | No stall |   0.212|
| YJL098W   | SAP185    | No stall |   0.211|
| YPL115C   | BEM3      | Stall    |   0.211|
| YKR020W   | VPS51     | No stall |   0.211|
| YMR233W   | TRI1      | Stall    |   0.211|
| YOL108C   | INO4      | Stall    |   0.211|
| YJR127C   | RSF2      | Stall    |   0.210|
| YHR183W   | GND1      | No stall |   0.209|
| YDR007W   | TRP1      | No stall |   0.209|
| YDR500C   | RPL37B    | No stall |   0.209|
| YKL186C   | MTR2      | No stall |   0.209|
| YOR159C   | SME1      | No stall |   0.209|
| YCL050C   | APA1      | No stall |   0.208|
| YOR370C   | MRS6      | No stall |   0.208|
| YDR349C   | YPS7      | No stall |   0.208|
| YDL105W   | NSE4      | Stall    |   0.207|
| YLR270W   | DCS1      | No stall |   0.207|
| YOR034C   | AKR2      | No stall |   0.207|
| YPL051W   | ARL3      | No stall |   0.206|
| YPL232W   | SSO1      | Stall    |   0.206|
| YML032C   | RAD52     | No stall |   0.205|
| YDL207W   | GLE1      | Stall    |   0.205|
| YHL015W   | RPS20     | No stall |   0.204|
| YBR255W   | MTC4      | Stall    |   0.204|
| YNL020C   | ARK1      | No stall |   0.204|
| YMR309C   | NIP1      | No stall |   0.203|
| YJL076W   | NET1      | Stall    |   0.203|
| YIL118W   | RHO3      | No stall |   0.203|
| YDR280W   | RRP45     | No stall |   0.201|
| YKL015W   | PUT3      | No stall |   0.201|
| YDR234W   | LYS4      | No stall |   0.201|
| YJL184W   | GON7      | No stall |   0.201|
| YKR093W   | PTR2      | No stall |   0.201|
| YLL034C   | RIX7      | Stall    |   0.200|
| YBR136W   | MEC1      | No stall |   0.200|
| YDL184C   | RPL41A    | Stall    |   0.200|
| YBR162C   | TOS1      | No stall |   0.200|
| YKL058W   | TOA2      | No stall |   0.199|
| YBL024W   | NCL1      | Stall    |   0.199|
| YBR084W   | MIS1      | No stall |   0.199|
| YJL004C   | SYS1      | No stall |   0.198|
| YDR302W   | GPI11     | Stall    |   0.198|
| YBR297W   | MAL33     | No stall |   0.197|
| YDR459C   | PFA5      | Stall    |   0.197|
| YDR060W   | MAK21     | Stall    |   0.197|
| YDL140C   | RPO21     | Stall    |   0.197|
| YOR275C   | RIM20     | No stall |   0.196|
| YJL129C   | TRK1      | Stall    |   0.196|
| YHR108W   | GGA2      | No stall |   0.196|
| YKL204W   | EAP1      | Stall    |   0.196|
| YPR031W   | NTO1      | No stall |   0.195|
| YLR380W   | CSR1      | No stall |   0.195|
| YHR123W   | EPT1      | No stall |   0.194|
| YNL191W   | DUG3      | No stall |   0.194|
| YPL101W   | ELP4      | No stall |   0.194|
| YDR487C   | RIB3      | No stall |   0.194|
| YDL025C   | RTK1      | Stall    |   0.193|
| YOR242C   | SSP2      | Stall    |   0.193|
| YNL240C   | NAR1      | No stall |   0.193|
| YBL021C   | HAP3      | No stall |   0.193|
| YPR080W   | TEF1      | No stall |   0.193|
| YDL148C   | NOP14     | Stall    |   0.193|
| YOL122C   | SMF1      | No stall |   0.192|
| YKL150W   | MCR1      | No stall |   0.192|
| YKR074W   | AIM29     | No stall |   0.191|
| YML023C   | NSE5      | No stall |   0.191|
| YGL223C   | COG1      | No stall |   0.191|
| YLR398C   | SKI2      | No stall |   0.191|
| YCL051W   | LRE1      | Stall    |   0.191|
| YJR044C   | VPS55     | No stall |   0.191|
| YEL059C-A | SOM1      | No stall |   0.190|
| YBL055C   | NA        | No stall |   0.189|
| YPL269W   | KAR9      | No stall |   0.189|
| YBL088C   | TEL1      | No stall |   0.188|
| YFR034C   | PHO4      | No stall |   0.188|
| YPR156C   | TPO3      | Stall    |   0.188|
| YIL155C   | GUT2      | No stall |   0.187|
| YDR389W   | SAC7      | Stall    |   0.187|
| YNL268W   | LYP1      | No stall |   0.187|
| YIL051C   | MMF1      | No stall |   0.187|
| YBR104W   | YMC2      | No stall |   0.186|
| YKR048C   | NAP1      | No stall |   0.186|
| YPL075W   | GCR1      | No stall |   0.185|
| YIL115C   | NUP159    | No stall |   0.185|
| YEL001C   | IRC22     | No stall |   0.185|
| YGR154C   | GTO1      | No stall |   0.185|
| YNL099C   | OCA1      | No stall |   0.185|
| YFR029W   | PTR3      | Stall    |   0.184|
| YGR193C   | PDX1      | No stall |   0.184|
| YPR162C   | ORC4      | No stall |   0.184|
| YOR310C   | NOP58     | Stall    |   0.184|
| YDR151C   | CTH1      | Stall    |   0.183|
| YDR147W   | EKI1      | No stall |   0.183|
| YAL023C   | PMT2      | No stall |   0.183|
| YPL198W   | RPL7B     | No stall |   0.183|
| YDR023W   | SES1      | No stall |   0.182|
| YLR098C   | CHA4      | Stall    |   0.182|
| YOR176W   | HEM15     | No stall |   0.182|
| YGL208W   | SIP2      | No stall |   0.181|
| YGR239C   | PEX21     | No stall |   0.181|
| YHR067W   | HTD2      | No stall |   0.181|
| YGR284C   | ERV29     | No stall |   0.181|
| YEL016C   | NPP2      | No stall |   0.181|
| YGL013C   | PDR1      | Stall    |   0.179|
| YER025W   | GCD11     | No stall |   0.179|
| YBR092C   | PHO3      | No stall |   0.179|
| YLR360W   | VPS38     | No stall |   0.179|
| YGR120C   | COG2      | No stall |   0.179|
| YOR260W   | GCD1      | No stall |   0.179|
| YJL081C   | ARP4      | No stall |   0.178|
| YLR353W   | BUD8      | No stall |   0.178|
| YIL114C   | POR2      | No stall |   0.178|
| YOR116C   | RPO31     | No stall |   0.178|
| YEL006W   | YEA6      | No stall |   0.178|
| YPL070W   | MUK1      | No stall |   0.177|
| YFR019W   | FAB1      | Stall    |   0.177|
| YHR137W   | ARO9      | No stall |   0.176|
| YNL136W   | EAF7      | Stall    |   0.176|
| YDR312W   | SSF2      | Stall    |   0.176|
| YJL093C   | TOK1      | No stall |   0.175|
| YNL322C   | KRE1      | No stall |   0.175|
| YKR090W   | PXL1      | No stall |   0.175|
| YAR014C   | BUD14     | No stall |   0.175|
| YML111W   | BUL2      | Stall    |   0.175|
| YLL003W   | SFI1      | No stall |   0.174|
| YJR112W   | NNF1      | No stall |   0.174|
| YDR464W   | SPP41     | Stall    |   0.174|
| YGR072W   | UPF3      | Stall    |   0.173|
| YCR046C   | IMG1      | Stall    |   0.173|
| YNR039C   | ZRG17     | Stall    |   0.173|
| YOR311C   | DGK1      | No stall |   0.173|
| YJL068C   | NA        | No stall |   0.173|
| YGR250C   | RIE1      | No stall |   0.173|
| YKR056W   | TRM2      | Stall    |   0.172|
| YDR038C   | ENA5      | No stall |   0.171|
| YDR545W   | YRF1-5    | No stall |   0.171|
| YPR058W   | YMC1      | No stall |   0.171|
| YFL049W   | SWP82     | Stall    |   0.171|
| YPL270W   | MDL2      | No stall |   0.171|
| YIL095W   | PRK1      | Stall    |   0.170|
| YOL070C   | NBA1      | Stall    |   0.170|
| YCL059C   | KRR1      | Stall    |   0.170|
| YIL083C   | CAB2      | No stall |   0.170|
| YLR095C   | IOC2      | Stall    |   0.170|
| YAL044C   | GCV3      | No stall |   0.169|
| YPR169W   | JIP5      | Stall    |   0.169|
| YDR419W   | RAD30     | No stall |   0.169|
| YLR068W   | FYV7      | Stall    |   0.169|
| YPL175W   | SPT14     | No stall |   0.168|
| YKR037C   | SPC34     | No stall |   0.168|
| YKR103W   | NFT1      | No stall |   0.168|
| YMR055C   | BUB2      | No stall |   0.168|
| YER056C-A | RPL34A    | Stall    |   0.168|
| YFL034C-A | RPL22B    | No stall |   0.167|
| YKL212W   | SAC1      | No stall |   0.167|
| YDR039C   | ENA2      | No stall |   0.167|
| YOL003C   | PFA4      | No stall |   0.166|
| YOR129C   | AFI1      | No stall |   0.166|
| YDR399W   | HPT1      | No stall |   0.166|
| YCR052W   | RSC6      | No stall |   0.165|
| YMR301C   | ATM1      | No stall |   0.165|
| YJL062W-A | COA3      | No stall |   0.164|
| YOR048C   | RAT1      | Stall    |   0.164|
| YGL147C   | RPL9A     | No stall |   0.164|
| YGL091C   | NBP35     | No stall |   0.164|
| YCL043C   | PDI1      | No stall |   0.164|
| YGR296W   | YRF1-7    | No stall |   0.163|
| YDL112W   | TRM3      | No stall |   0.163|
| YOR085W   | OST3      | No stall |   0.162|
| YGR119C   | NUP57     | No stall |   0.162|
| YMR278W   | PRM15     | No stall |   0.162|
| YPR127W   | NA        | No stall |   0.162|
| YDL133C-A | RPL41A    | Stall    |   0.162|
| YPL174C   | NIP100    | Stall    |   0.161|
| YJL106W   | IME2      | No stall |   0.161|
| YOL095C   | HMI1      | No stall |   0.161|
| YHR143W-A | RPC10     | No stall |   0.161|
| YPL133C   | RDS2      | No stall |   0.160|
| YLR467W   | YRF1-5    | No stall |   0.160|
| YDR047W   | HEM12     | No stall |   0.160|
| YNL087W   | TCB2      | No stall |   0.160|
| YBR160W   | CDC28     | No stall |   0.160|
| YMR314W   | PRE5      | No stall |   0.160|
| YGL008C   | PMA1      | No stall |   0.159|
| YHR077C   | NMD2      | No stall |   0.159|
| YPL047W   | SGF11     | No stall |   0.159|
| YNL138W   | SRV2      | Stall    |   0.158|
| YER029C   | SMB1      | Stall    |   0.158|
| YJL110C   | GZF3      | No stall |   0.158|
| YDR037W   | KRS1      | No stall |   0.158|
| YNR052C   | POP2      | No stall |   0.158|
| YIL050W   | PCL7      | No stall |   0.158|
| YKL144C   | RPC25     | No stall |   0.158|
| YMR263W   | SAP30     | Stall    |   0.158|
| YPL005W   | AEP3      | Stall    |   0.158|
| YDR065W   | RRG1      | No stall |   0.157|
| YOL021C   | DIS3      | No stall |   0.157|
| YDR408C   | ADE8      | No stall |   0.157|
| YHL013C   | OTU2      | Stall    |   0.157|
| YGL198W   | YIP4      | No stall |   0.156|
| YMR189W   | GCV2      | No stall |   0.156|
| YBL035C   | POL12     | No stall |   0.156|
| YEL004W   | YEA4      | No stall |   0.156|
| YDR069C   | DOA4      | No stall |   0.154|
| YKR079C   | TRZ1      | No stall |   0.154|
| YHL020C   | OPI1      | No stall |   0.154|
| YER172C   | BRR2      | No stall |   0.154|
| YER139C   | RTR1      | No stall |   0.153|
| YNL264C   | PDR17     | No stall |   0.153|
| YMR073C   | IRC21     | No stall |   0.153|
| YPR186C   | PZF1      | Stall    |   0.153|
| YJL198W   | PHO90     | No stall |   0.152|
| YOR362C   | PRE10     | No stall |   0.152|
| YIL122W   | POG1      | No stall |   0.151|
| YLR320W   | MMS22     | Stall    |   0.151|
| YMR119W   | ASI1      | No stall |   0.151|
| YHR098C   | SFB3      | No stall |   0.151|
| YAR003W   | SWD1      | No stall |   0.150|
| YOR047C   | STD1      | No stall |   0.150|
| YMR021C   | MAC1      | No stall |   0.150|
| YIL056W   | VHR1      | Stall    |   0.150|
| YER063W   | THO1      | No stall |   0.149|
| YML013W   | UBX2      | No stall |   0.149|
| YNL002C   | RLP7      | Stall    |   0.149|
| YDR239C   | NA        | No stall |   0.148|
| YML098W   | TAF13     | No stall |   0.147|
| YKL205W   | LOS1      | No stall |   0.147|
| YNL253W   | TEX1      | No stall |   0.147|
| YHR003C   | TCD1      | Stall    |   0.147|
| YLR234W   | TOP3      | No stall |   0.147|
| YJL125C   | GCD14     | No stall |   0.147|
| YMR313C   | TGL3      | No stall |   0.146|
| YHR186C   | KOG1      | No stall |   0.146|
| YCL005W-A | VMA9      | No stall |   0.146|
| YGR257C   | MTM1      | No stall |   0.146|
| YAL051W   | OAF1      | No stall |   0.146|
| YJL178C   | ATG27     | No stall |   0.145|
| YML007W   | YAP1      | Stall    |   0.145|
| YJL061W   | NUP82     | No stall |   0.145|
| YHR103W   | SBE22     | No stall |   0.145|
| YJL094C   | KHA1      | No stall |   0.145|
| YDR178W   | SDH4      | No stall |   0.144|
| YKR009C   | FOX2      | No stall |   0.144|
| YGR003W   | CUL3      | No stall |   0.144|
| YCR073C   | SSK22     | No stall |   0.144|
| YDR052C   | DBF4      | Stall    |   0.144|
| YMR127C   | SAS2      | No stall |   0.143|
| YLR385C   | SWC7      | No stall |   0.143|
| YHR079C   | IRE1      | Stall    |   0.143|
| YER089C   | PTC2      | No stall |   0.142|
| YBL106C   | SRO77     | No stall |   0.142|
| YMR120C   | ADE17     | No stall |   0.141|
| YDL139C   | SCM3      | Stall    |   0.141|
| YKR036C   | CAF4      | No stall |   0.141|
| YNL183C   | NPR1      | No stall |   0.140|
| YHR082C   | KSP1      | Stall    |   0.140|
| YER021W   | RPN3      | No stall |   0.140|
| YER070W   | RNR1      | No stall |   0.140|
| YML064C   | TEM1      | No stall |   0.140|
| YPL147W   | PXA1      | Stall    |   0.139|
| YHR065C   | RRP3      | Stall    |   0.139|
| YOL001W   | PHO80     | No stall |   0.139|
| YBR118W   | TEF1      | No stall |   0.138|
| YOL022C   | TSR4      | No stall |   0.138|
| YGR041W   | BUD9      | No stall |   0.138|
| YLR405W   | DUS4      | No stall |   0.138|
| YJR090C   | GRR1      | No stall |   0.138|
| YLL031C   | GPI13     | No stall |   0.137|
| YER078C   | ICP55     | No stall |   0.137|
| YBL052C   | SAS3      | Stall    |   0.137|
| YNL321W   | VNX1      | No stall |   0.137|
| YER007C-A | TMA20     | No stall |   0.137|
| YHR170W   | NMD3      | Stall    |   0.137|
| YLR022C   | SDO1      | Stall    |   0.137|
| YJR010W   | MET3      | No stall |   0.136|
| YOL063C   | CRT10     | No stall |   0.136|
| YIL010W   | DOT5      | No stall |   0.136|
| YOR330C   | MIP1      | Stall    |   0.136|
| YDR332W   | IRC3      | No stall |   0.135|
| YLR244C   | MAP1      | Stall    |   0.135|
| YHL032C   | GUT1      | No stall |   0.134|
| YKR030W   | GMH1      | No stall |   0.134|
| YDL054C   | MCH1      | No stall |   0.134|
| YDR321W   | ASP1      | No stall |   0.134|
| YOR138C   | RUP1      | Stall    |   0.133|
| YDR019C   | GCV1      | No stall |   0.133|
| YNL061W   | NOP2      | Stall    |   0.132|
| YDR486C   | VPS60     | No stall |   0.132|
| YOR241W   | MET7      | No stall |   0.132|
| YFR040W   | SAP155    | Stall    |   0.132|
| YLR048W   | RPS0B     | No stall |   0.132|
| YBR153W   | RIB7      | No stall |   0.132|
| YER173W   | RAD24     | Stall    |   0.132|
| YPL152W   | RRD2      | No stall |   0.132|
| YOR295W   | UAF30     | No stall |   0.132|
| YNR024W   | MPP6      | Stall    |   0.132|
| YNL024C-A | KSH1      | No stall |   0.131|
| YNL159C   | ASI2      | No stall |   0.131|
| YGR089W   | NNF2      | Stall    |   0.130|
| YBR005W   | RCR1      | Stall    |   0.130|
| YHR085W   | IPI1      | Stall    |   0.130|
| YMR214W   | SCJ1      | No stall |   0.129|
| YOR127W   | RGA1      | Stall    |   0.128|
| YML124C   | TUB3      | No stall |   0.128|
| YNL106C   | INP52     | No stall |   0.128|
| YOR335C   | ALA1      | No stall |   0.127|
| YAL059W   | ECM1      | Stall    |   0.127|
| YLR016C   | PML1      | No stall |   0.127|
| YEL034W   | HYP2      | No stall |   0.127|
| YPL226W   | NEW1      | Stall    |   0.126|
| YIL162W   | SUC2      | No stall |   0.126|
| YBR243C   | ALG7      | No stall |   0.125|
| YMR307W   | GAS1      | No stall |   0.125|
| YHR188C   | GPI16     | No stall |   0.125|
| YDL100C   | GET3      | No stall |   0.125|
| YIR018W   | YAP5      | No stall |   0.125|
| YDL064W   | UBC9      | No stall |   0.125|
| YDR176W   | NGG1      | Stall    |   0.124|
| YOR092W   | ECM3      | No stall |   0.124|
| YFR038W   | IRC5      | Stall    |   0.124|
| YOR126C   | IAH1      | No stall |   0.123|
| YJR084W   | NA        | Stall    |   0.123|
| YHR181W   | SVP26     | No stall |   0.123|
| YGR199W   | PMT6      | No stall |   0.123|
| YOR144C   | ELG1      | No stall |   0.122|
| YOR182C   | RPS30A    | Stall    |   0.122|
| YOR213C   | SAS5      | No stall |   0.122|
| YJR076C   | CDC11     | Stall    |   0.122|
| YDR488C   | PAC11     | No stall |   0.122|
| YDR026C   | NSI1      | Stall    |   0.121|
| YEL071W   | DLD3      | No stall |   0.121|
| YCR020W-B | HTL1      | No stall |   0.120|
| YIL129C   | TAO3      | No stall |   0.120|
| YBL091C   | MAP2      | Stall    |   0.120|
| YER131W   | RPS26B    | No stall |   0.120|
| YPL194W   | DDC1      | No stall |   0.119|
| YAL005C   | SSA1      | No stall |   0.119|
| YPL221W   | FLC1      | No stall |   0.119|
| YCR068W   | ATG15     | No stall |   0.118|
| YNL148C   | ALF1      | No stall |   0.118|
| YHR146W   | CRP1      | Stall    |   0.118|
| YOR209C   | NPT1      | No stall |   0.118|
| YGL075C   | MPS2      | No stall |   0.118|
| YDR323C   | PEP7      | No stall |   0.117|
| YMR275C   | BUL1      | Stall    |   0.117|
| YHR194W   | MDM31     | No stall |   0.117|
| YDR148C   | KGD2      | No stall |   0.117|
| YPL256C   | CLN2      | No stall |   0.116|
| YJR068W   | RFC2      | No stall |   0.116|
| YGR204W   | ADE3      | No stall |   0.116|
| YDL127W   | PCL2      | No stall |   0.116|
| YNL255C   | GIS2      | No stall |   0.116|
| YOL027C   | MDM38     | No stall |   0.116|
| YEL065W   | SIT1      | No stall |   0.115|
| YMR208W   | ERG12     | No stall |   0.115|
| YDR449C   | UTP6      | Stall    |   0.115|
| YOR107W   | RGS2      | No stall |   0.114|
| YLR064W   | PER33     | No stall |   0.114|
| YJR052W   | RAD7      | Stall    |   0.114|
| YJL128C   | PBS2      | No stall |   0.114|
| YKL110C   | KTI12     | No stall |   0.114|
| YGR229C   | SMI1      | No stall |   0.114|
| YDL141W   | BPL1      | No stall |   0.113|
| YPR024W   | YME1      | No stall |   0.113|
| YNL064C   | YDJ1      | No stall |   0.113|
| YER099C   | PRS2      | No stall |   0.112|
| YIL084C   | SDS3      | Stall    |   0.112|
| YJL111W   | CCT7      | No stall |   0.112|
| YER147C   | SCC4      | No stall |   0.112|
| YAL015C   | NTG1      | No stall |   0.111|
| YPR125W   | YLH47     | No stall |   0.111|
| YNL209W   | SSB2      | No stall |   0.110|
| YER028C   | MIG3      | Stall    |   0.110|
| YBR065C   | ECM2      | Stall    |   0.110|
| YOR217W   | RFC1      | Stall    |   0.110|
| YGL131C   | SNT2      | Stall    |   0.110|
| YPR185W   | ATG13     | No stall |   0.110|
| YDL031W   | DBP10     | Stall    |   0.109|
| YKL138C-A | HSK3      | No stall |   0.109|
| YDR184C   | ATC1      | No stall |   0.109|
| YMR128W   | ECM16     | Stall    |   0.109|
| YLR427W   | MAG2      | Stall    |   0.109|
| YPL151C   | PRP46     | No stall |   0.108|
| YFR004W   | RPN11     | No stall |   0.108|
| YAL041W   | CDC24     | No stall |   0.108|
| YDL077C   | VAM6      | No stall |   0.108|
| YER049W   | TPA1      | No stall |   0.107|
| YJR064W   | CCT5      | No stall |   0.107|
| YOR051C   | ETT1      | Stall    |   0.107|
| YER003C   | PMI40     | No stall |   0.107|
| YJL207C   | LAA1      | Stall    |   0.106|
| YFL031W   | HAC1      | Stall    |   0.106|
| YLR075W   | RPL10     | No stall |   0.106|
| YJL222W   | VTH2      | Stall    |   0.106|
| YBR015C   | MNN2      | Stall    |   0.105|
| YGL145W   | TIP20     | No stall |   0.105|
| YBL004W   | UTP20     | Stall    |   0.105|
| YOL066C   | RIB2      | No stall |   0.105|
| YKL050C   | NA        | No stall |   0.104|
| YGL106W   | MLC1      | No stall |   0.104|
| YFR013W   | IOC3      | Stall    |   0.104|
| YDR320C   | SWA2      | No stall |   0.104|
| YLR430W   | SEN1      | Stall    |   0.104|
| YGR143W   | SKN1      | No stall |   0.103|
| YPL119C   | DBP1      | No stall |   0.103|
| YGL105W   | ARC1      | No stall |   0.103|
| YPL030W   | TRM44     | No stall |   0.103|
| YPR067W   | ISA2      | No stall |   0.102|
| YDR081C   | PDC2      | Stall    |   0.102|
| YPL103C   | FMP30     | No stall |   0.102|
| YFL010C   | WWM1      | No stall |   0.102|
| YLR174W   | IDP2      | No stall |   0.101|
| YHR091C   | MSR1      | No stall |   0.101|
| YML048W   | GSF2      | No stall |   0.100|
| YLR393W   | ATP10     | No stall |   0.100|
| YDL175C   | AIR2      | No stall |   0.100|
| YLR410W   | VIP1      | Stall    |   0.100|
| YKR082W   | NUP133    | No stall |   0.100|
| YGR175C   | ERG1      | No stall |   0.100|
| YER190W   | YRF1-2    | No stall |   0.099|
| YER098W   | UBP9      | Stall    |   0.099|
| YHR114W   | BZZ1      | No stall |   0.098|
| YOR291W   | YPK9      | Stall    |   0.098|
| YJL091C   | GWT1      | No stall |   0.098|
| YMR129W   | POM152    | No stall |   0.098|
| YML080W   | DUS1      | No stall |   0.097|
| YER061C   | CEM1      | No stall |   0.097|
| YJL189W   | RPL39     | Stall    |   0.097|
| YDL188C   | PPH22     | No stall |   0.097|
| YMR227C   | TAF7      | Stall    |   0.097|
| YPR048W   | TAH18     | No stall |   0.097|
| YOL089C   | HAL9      | Stall    |   0.096|
| YLR404W   | SEI1      | No stall |   0.096|
| YFL001W   | DEG1      | Stall    |   0.096|
| YGL232W   | TAN1      | No stall |   0.095|
| YDR421W   | ARO80     | No stall |   0.095|
| YDL198C   | GGC1      | No stall |   0.095|
| YJR009C   | TDH2      | No stall |   0.094|
| YOR001W   | RRP6      | Stall    |   0.094|
| YER092W   | IES5      | No stall |   0.094|
| YNL172W   | APC1      | No stall |   0.094|
| YOR160W   | MTR10     | No stall |   0.093|
| YEL017W   | GTT3      | No stall |   0.093|
| YLR299W   | ECM38     | No stall |   0.093|
| YBR110W   | ALG1      | No stall |   0.092|
| YMR223W   | UBP8      | No stall |   0.092|
| YGR009C   | SEC9      | Stall    |   0.092|
| YOR149C   | SMP3      | No stall |   0.092|
| YOR194C   | TOA1      | No stall |   0.092|
| YJL145W   | SFH5      | No stall |   0.091|
| YPL149W   | ATG5      | No stall |   0.091|
| YDR353W   | TRR1      | No stall |   0.091|
| YMR232W   | FUS2      | No stall |   0.091|
| YHR039C-A | VMA10     | No stall |   0.091|
| YPL046C   | ELC1      | No stall |   0.091|
| YJR060W   | CBF1      | No stall |   0.091|
| YBL079W   | NUP170    | No stall |   0.091|
| YOR201C   | MRM1      | No stall |   0.091|
| YGL111W   | NSA1      | Stall    |   0.091|
| YKR010C   | TOF2      | Stall    |   0.091|
| YPR069C   | SPE3      | No stall |   0.090|
| YOL080C   | REX4      | No stall |   0.090|
| YOR336W   | KRE5      | No stall |   0.090|
| YML057W   | CMP2      | No stall |   0.090|
| YDL143W   | CCT4      | No stall |   0.090|
| YIR015W   | RPR2      | No stall |   0.088|
| YFL009W   | CDC4      | No stall |   0.088|
| YLR035C   | MLH2      | No stall |   0.088|
| YOR042W   | CUE5      | Stall    |   0.088|
| YIL006W   | YIA6      | No stall |   0.088|
| YAR031W   | PRM9      | No stall |   0.088|
| YDL115C   | IWR1      | No stall |   0.087|
| YPL243W   | SRP68     | Stall    |   0.087|
| YLR305C   | STT4      | Stall    |   0.087|
| YNR051C   | BRE5      | No stall |   0.087|
| YIL153W   | RRD1      | No stall |   0.087|
| YBL101C   | ECM21     | No stall |   0.087|
| YIL128W   | MET18     | No stall |   0.086|
| YGR281W   | YOR1      | Stall    |   0.086|
| YGR238C   | KEL2      | No stall |   0.086|
| YBL097W   | BRN1      | Stall    |   0.086|
| YIR009W   | MSL1      | No stall |   0.086|
| YGL219C   | MDM34     | Stall    |   0.085|
| YNL271C   | BNI1      | Stall    |   0.085|
| YLR466W   | YRF1-4    | No stall |   0.084|
| YER036C   | ARB1      | Stall    |   0.084|
| YDR083W   | RRP8      | Stall    |   0.084|
| YDR117C   | TMA64     | Stall    |   0.084|
| YLR024C   | UBR2      | Stall    |   0.083|
| YJL073W   | JEM1      | No stall |   0.083|
| YDL205C   | HEM3      | No stall |   0.083|
| YCR035C   | RRP43     | No stall |   0.083|
| YLR363C   | NMD4      | No stall |   0.082|
| YNR008W   | LRO1      | Stall    |   0.082|
| YDR347W   | MRP1      | No stall |   0.082|
| YEL054C   | RPL12B    | Stall    |   0.082|
| YJL203W   | PRP21     | Stall    |   0.082|
| YDL049C   | KNH1      | No stall |   0.081|
| YPL059W   | GRX5      | No stall |   0.081|
| YIL173W   | VTH1      | Stall    |   0.080|
| YDR346C   | SVF1      | No stall |   0.080|
| YPR145W   | ASN1      | No stall |   0.080|
| YGL003C   | CDH1      | No stall |   0.080|
| YBL089W   | AVT5      | No stall |   0.079|
| YIL079C   | AIR1      | No stall |   0.079|
| YML024W   | RPS17A    | Stall    |   0.079|
| YOR207C   | RET1      | No stall |   0.079|
| YKL046C   | DCW1      | No stall |   0.078|
| YPR095C   | SYT1      | Stall    |   0.078|
| YJL208C   | NUC1      | No stall |   0.078|
| YBR114W   | RAD16     | Stall    |   0.077|
| YOR353C   | SOG2      | No stall |   0.077|
| YDL070W   | BDF2      | Stall    |   0.077|
| YGL086W   | MAD1      | No stall |   0.077|
| YMR156C   | TPP1      | No stall |   0.076|
| YPR180W   | AOS1      | No stall |   0.076|
| YMR216C   | SKY1      | Stall    |   0.076|
| YBR011C   | IPP1      | No stall |   0.076|
| YOR344C   | TYE7      | Stall    |   0.075|
| YBR170C   | NPL4      | No stall |   0.075|
| YPL023C   | MET12     | No stall |   0.075|
| YDR400W   | URH1      | No stall |   0.075|
| YCL045C   | EMC1      | No stall |   0.075|
| YER175C   | TMT1      | No stall |   0.075|
| YPR086W   | SUA7      | No stall |   0.073|
| YHR080C   | LAM4      | Stall    |   0.073|
| YGL110C   | CUE3      | Stall    |   0.073|
| YOL062C   | APM4      | No stall |   0.073|
| YLR190W   | MMR1      | No stall |   0.073|
| YCR054C   | CTR86     | No stall |   0.073|
| YDR001C   | NTH1      | No stall |   0.072|
| YGL150C   | INO80     | Stall    |   0.072|
| YKL197C   | PEX1      | No stall |   0.072|
| YFR037C   | RSC8      | No stall |   0.072|
| YBL105C   | PKC1      | Stall    |   0.072|
| YJL217W   | REE1      | No stall |   0.071|
| YKR007W   | MEH1      | No stall |   0.071|
| YGR070W   | ROM1      | No stall |   0.071|
| YER168C   | CCA1      | No stall |   0.071|
| YHR133C   | NSG1      | No stall |   0.070|
| YPL125W   | KAP120    | No stall |   0.070|
| YNR011C   | PRP2      | No stall |   0.070|
| YGL216W   | KIP3      | No stall |   0.069|
| YNR055C   | HOL1      | No stall |   0.069|
| YBL042C   | FUI1      | No stall |   0.069|
| YGL029W   | CGR1      | Stall    |   0.069|
| YPL172C   | COX10     | No stall |   0.069|
| YML056C   | IMD4      | No stall |   0.069|
| YDL047W   | SIT4      | No stall |   0.068|
| YJL026W   | RNR2      | No stall |   0.068|
| YOR297C   | TIM18     | No stall |   0.068|
| YCR053W   | THR4      | No stall |   0.068|
| YPL233W   | NSL1      | No stall |   0.067|
| YDL015C   | TSC13     | Stall    |   0.067|
| YIL088C   | AVT7      | No stall |   0.066|
| YLR113W   | HOG1      | No stall |   0.066|
| YJR022W   | LSM8      | No stall |   0.066|
| YOR077W   | RTS2      | Stall    |   0.066|
| YHR012W   | VPS29     | No stall |   0.066|
| YDR368W   | YPR1      | No stall |   0.066|
| YGL011C   | SCL1      | No stall |   0.066|
| YGR270W   | YTA7      | Stall    |   0.065|
| YNL003C   | PET8      | No stall |   0.065|
| YDR135C   | YCF1      | No stall |   0.064|
| YLR019W   | PSR2      | No stall |   0.063|
| YGL186C   | TPN1      | Stall    |   0.063|
| YLR246W   | ERF2      | No stall |   0.062|
| YOR078W   | BUD21     | Stall    |   0.062|
| YML125C   | PGA3      | No stall |   0.062|
| YJL097W   | PHS1      | No stall |   0.062|
| YLR141W   | RRN5      | No stall |   0.062|
| YOL045W   | PSK2      | No stall |   0.062|
| YMR032W   | HOF1      | No stall |   0.062|
| YLR175W   | CBF5      | Stall    |   0.061|
| YGR188C   | BUB1      | Stall    |   0.061|
| YDR129C   | SAC6      | No stall |   0.061|
| YLL021W   | SPA2      | No stall |   0.061|
| YHR118C   | ORC6      | No stall |   0.061|
| YNL045W   | LAP2      | No stall |   0.060|
| YMR238W   | DFG5      | No stall |   0.060|
| YNL229C   | URE2      | No stall |   0.060|
| YGL012W   | ERG4      | No stall |   0.060|
| YDR381W   | YRA1      | Stall    |   0.060|
| YKL057C   | NUP120    | No stall |   0.060|
| YMR176W   | ECM5      | No stall |   0.060|
| YKL129C   | MYO3      | Stall    |   0.059|
| YDR385W   | EFT1      | No stall |   0.059|
| YLR102C   | APC9      | No stall |   0.058|
| YOR184W   | SER1      | No stall |   0.058|
| YBR272C   | HSM3      | No stall |   0.058|
| YIL064W   | EFM4      | No stall |   0.058|
| YHR165C   | PRP8      | Stall    |   0.058|
| YLR005W   | SSL1      | Stall    |   0.058|
| YGR241C   | YAP1802   | No stall |   0.058|
| YNL316C   | PHA2      | No stall |   0.058|
| YNL091W   | NST1      | Stall    |   0.057|
| YLL032C   | NA        | No stall |   0.057|
| YDR466W   | PKH3      | Stall    |   0.057|
| YLR249W   | YEF3      | Stall    |   0.056|
| YKR086W   | PRP16     | Stall    |   0.056|
| YML025C   | YML6      | No stall |   0.056|
| YGR004W   | PEX31     | Stall    |   0.056|
| YGL148W   | ARO2      | No stall |   0.056|
| YDL020C   | RPN4      | No stall |   0.056|
| YFL048C   | EMP47     | No stall |   0.055|
| YOR094W   | ARF3      | No stall |   0.054|
| YDR162C   | NBP2      | No stall |   0.054|
| YKL139W   | CTK1      | No stall |   0.054|
| YJL095W   | BCK1      | No stall |   0.053|
| YPL099C   | INA17     | No stall |   0.053|
| YJR045C   | SSC1      | No stall |   0.053|
| YNL328C   | MDJ2      | Stall    |   0.052|
| YJL148W   | RPA34     | Stall    |   0.052|
| YJL024C   | APS3      | No stall |   0.052|
| YPL227C   | ALG5      | No stall |   0.052|
| YGR192C   | TDH3      | No stall |   0.052|
| YNL139C   | THO2      | No stall |   0.052|
| YAL029C   | MYO4      | Stall    |   0.051|
| YGR166W   | TRS65     | No stall |   0.051|
| YHR073W   | OSH3      | Stall    |   0.051|
| YJR042W   | NUP85     | No stall |   0.051|
| YDR145W   | TAF12     | No stall |   0.051|
| YPL028W   | ERG10     | No stall |   0.051|
| YLR038C   | COX12     | No stall |   0.051|
| YHR132C   | ECM14     | No stall |   0.050|
| YHR158C   | KEL1      | No stall |   0.050|
| YMR036C   | MIH1      | No stall |   0.050|
| YCL014W   | BUD3      | No stall |   0.049|
| YDL076C   | RXT3      | No stall |   0.049|
| YGL017W   | ATE1      | No stall |   0.049|
| YDR341C   | NA        | No stall |   0.049|
| YPL249C   | GYP5      | No stall |   0.049|
| YPL262W   | FUM1      | No stall |   0.048|
| YJR123W   | RPS5      | No stall |   0.048|
| YFL062W   | COS4      | No stall |   0.048|
| YJL101C   | GSH1      | No stall |   0.048|
| YDR418W   | RPL12B    | Stall    |   0.048|
| YJL042W   | MHP1      | No stall |   0.048|
| YDR198C   | RKM2      | No stall |   0.048|
| YLR295C   | ATP14     | No stall |   0.047|
| YDL108W   | KIN28     | No stall |   0.047|
| YGL151W   | NUT1      | No stall |   0.047|
| YMR246W   | FAA4      | No stall |   0.047|
| YMR272C   | SCS7      | No stall |   0.047|
| YOL025W   | LAG2      | No stall |   0.046|
| YKR031C   | SPO14     | Stall    |   0.046|
| YAL013W   | DEP1      | No stall |   0.045|
| YPL212C   | PUS1      | Stall    |   0.045|
| YNR074C   | AIF1      | No stall |   0.045|
| YFL007W   | BLM10     | Stall    |   0.045|
| YDR104C   | SPO71     | No stall |   0.045|
| YOR019W   | NA        | Stall    |   0.045|
| YNL201C   | PSY2      | No stall |   0.045|
| YNL278W   | CAF120    | No stall |   0.044|
| YDL007W   | RPT2      | Stall    |   0.044|
| YMR006C   | PLB2      | No stall |   0.044|
| YLL010C   | PSR1      | No stall |   0.043|
| YPL106C   | SSE1      | No stall |   0.043|
| YDL074C   | BRE1      | No stall |   0.043|
| YBR237W   | PRP5      | Stall    |   0.042|
| YPL146C   | NOP53     | Stall    |   0.042|
| YOR133W   | EFT1      | No stall |   0.041|
| YHR005C-A | TIM10     | No stall |   0.041|
| YCL040W   | GLK1      | No stall |   0.041|
| YLR340W   | RPP0      | No stall |   0.041|
| YOR249C   | APC5      | No stall |   0.041|
| YKL108W   | SLD2      | Stall    |   0.040|
| YKL088W   | CAB3      | No stall |   0.040|
| YJL053W   | PEP8      | No stall |   0.040|
| YGL064C   | MRH4      | Stall    |   0.040|
| YPR036W   | VMA13     | No stall |   0.040|
| YPR189W   | SKI3      | No stall |   0.040|
| YCL052C   | PBN1      | Stall    |   0.040|
| YJL087C   | TRL1      | No stall |   0.040|
| YGR082W   | TOM20     | No stall |   0.040|
| YGL213C   | SKI8      | No stall |   0.040|
| YHL009C   | YAP3      | No stall |   0.039|
| YDR314C   | RAD34     | Stall    |   0.039|
| YNL039W   | BDP1      | Stall    |   0.039|
| YNL262W   | POL2      | Stall    |   0.039|
| YPL048W   | CAM1      | No stall |   0.038|
| YPL263C   | KEL3      | Stall    |   0.038|
| YJR049C   | UTR1      | No stall |   0.038|
| YPL110C   | GDE1      | No stall |   0.038|
| YER105C   | NUP157    | No stall |   0.038|
| YHR168W   | MTG2      | No stall |   0.037|
| YOR099W   | KTR1      | No stall |   0.037|
| YBR233W   | PBP2      | No stall |   0.036|
| YPL052W   | OAZ1      | No stall |   0.036|
| YGR002C   | SWC4      | No stall |   0.036|
| YDL090C   | RAM1      | No stall |   0.036|
| YNR047W   | FPK1      | Stall    |   0.036|
| YKR058W   | GLG1      | No stall |   0.036|
| YDR237W   | MRPL7     | No stall |   0.036|
| YPL098C   | MGR2      | No stall |   0.035|
| YBR060C   | ORC2      | Stall    |   0.035|
| YOR320C   | GNT1      | No stall |   0.035|
| YIL075C   | RPN2      | Stall    |   0.035|
| YGL127C   | SOH1      | No stall |   0.035|
| YEL015W   | EDC3      | No stall |   0.035|
| YBR130C   | SHE3      | No stall |   0.034|
| YPL020C   | ULP1      | Stall    |   0.033|
| YDR296W   | MHR1      | Stall    |   0.033|
| YIL048W   | NEO1      | No stall |   0.033|
| YER125W   | RSP5      | No stall |   0.033|
| YLL048C   | YBT1      | No stall |   0.033|
| YDL046W   | NPC2      | No stall |   0.033|
| YJL154C   | VPS35     | No stall |   0.032|
| YPL144W   | POC4      | No stall |   0.032|
| YJR040W   | GEF1      | No stall |   0.032|
| YML010W   | SPT5      | Stall    |   0.032|
| YOR067C   | ALG8      | No stall |   0.031|
| YAL027W   | SAW1      | No stall |   0.031|
| YKL209C   | STE6      | No stall |   0.031|
| YDR035W   | ARO3      | Stall    |   0.031|
| YKL175W   | ZRT3      | No stall |   0.030|
| YPL239W   | YAR1      | No stall |   0.030|
| YLR336C   | SGD1      | Stall    |   0.030|
| YOR308C   | SNU66     | Stall    |   0.030|
| YNR023W   | SNF12     | No stall |   0.030|
| YCR042C   | TAF2      | Stall    |   0.029|
| YLR242C   | ARV1      | No stall |   0.028|
| YPR075C   | OPY2      | Stall    |   0.028|
| YMR171C   | EAR1      | No stall |   0.028|
| YPR057W   | BRR1      | No stall |   0.028|
| YNR019W   | ARE2      | No stall |   0.027|
| YGL200C   | EMP24     | No stall |   0.027|
| YPR037C   | ERV2      | No stall |   0.027|
| YLR212C   | TUB4      | No stall |   0.027|
| YBR202W   | MCM7      | No stall |   0.027|
| YJR099W   | YUH1      | No stall |   0.026|
| YPL207W   | TYW1      | Stall    |   0.026|
| YGL202W   | ARO8      | No stall |   0.026|
| YER143W   | DDI1      | No stall |   0.025|
| YER157W   | COG3      | No stall |   0.025|
| YLR127C   | APC2      | No stall |   0.025|
| YGR063C   | SPT4      | No stall |   0.025|
| YGL119W   | COQ8      | No stall |   0.025|
| YKR028W   | SAP190    | No stall |   0.024|
| YHL014C   | YLF2      | No stall |   0.024|
| YCR026C   | NPP1      | No stall |   0.023|
| YDL236W   | PHO13     | No stall |   0.023|
| YER008C   | SEC3      | Stall    |   0.023|
| YBL007C   | SLA1      | Stall    |   0.023|
| YKL152C   | GPM1      | No stall |   0.023|
| YPR008W   | HAA1      | Stall    |   0.022|
| YBR169C   | SSE2      | No stall |   0.022|
| YDR079C-A | TFB5      | No stall |   0.022|
| YOL147C   | PEX11     | No stall |   0.022|
| YMR047C   | NUP116    | No stall |   0.021|
| YJL020C   | BBC1      | Stall    |   0.021|
| YKL059C   | MPE1      | No stall |   0.021|
| YGL060W   | YBP2      | No stall |   0.021|
| YHR047C   | AAP1      | No stall |   0.021|
| YBL019W   | APN2      | No stall |   0.021|
| YGR134W   | CAF130    | No stall |   0.020|
| YBR149W   | ARA1      | No stall |   0.020|
| YNL130C   | CPT1      | No stall |   0.020|
| YPL228W   | CET1      | No stall |   0.020|
| YJL096W   | MRPL49    | Stall    |   0.019|
| YJR100C   | AIM25     | No stall |   0.019|
| YHR195W   | NVJ1      | No stall |   0.019|
| YDL226C   | GCS1      | No stall |   0.019|
| YML015C   | TAF11     | No stall |   0.018|
| YCR066W   | RAD18     | Stall    |   0.018|
| YGL201C   | MCM6      | No stall |   0.017|
| YOR079C   | ATX2      | No stall |   0.017|
| YML086C   | ALO1      | Stall    |   0.017|
| YPL006W   | NCR1      | No stall |   0.016|
| YLR288C   | MEC3      | No stall |   0.016|
| YKR069W   | MET1      | No stall |   0.015|
| YLR085C   | ARP6      | No stall |   0.015|
| YML001W   | YPT7      | No stall |   0.015|
| YDL056W   | MBP1      | Stall    |   0.015|
| YBR044C   | TCM62     | No stall |   0.014|
| YNL125C   | ESBP6     | No stall |   0.014|
| YDR246W   | TRS23     | No stall |   0.014|
| YHR081W   | LRP1      | Stall    |   0.014|
| YMR218C   | TRS130    | Stall    |   0.014|
| YJR017C   | ESS1      | No stall |   0.014|
| YMR114C   | NA        | Stall    |   0.014|
| YGR221C   | TOS2      | No stall |   0.014|
| YHR164C   | DNA2      | Stall    |   0.013|
| YOR155C   | ISN1      | No stall |   0.013|
| YOR080W   | DIA2      | Stall    |   0.013|
| YPR199C   | ARR1      | Stall    |   0.013|
| YBR164C   | ARL1      | No stall |   0.013|
| YDL028C   | MPS1      | No stall |   0.012|
| YLR378C   | SEC61     | No stall |   0.012|
| YLR222C   | UTP13     | No stall |   0.012|
| YIL034C   | CAP2      | No stall |   0.011|
| YBL099W   | ATP1      | No stall |   0.011|
| YJR035W   | RAD26     | Stall    |   0.011|
| YNL326C   | PFA3      | No stall |   0.011|
| YDL165W   | CDC36     | No stall |   0.011|
| YCR093W   | CDC39     | No stall |   0.011|
| YKL009W   | MRT4      | No stall |   0.010|
| YDL229W   | SSB1      | No stall |   0.010|
| YNL065W   | AQR1      | No stall |   0.010|
| YNL049C   | SFB2      | No stall |   0.010|
| YPL141C   | FRK1      | Stall    |   0.009|
| YHR109W   | CTM1      | No stall |   0.009|
| YJR075W   | HOC1      | No stall |   0.009|
| YHR201C   | PPX1      | No stall |   0.009|
| YOR156C   | NFI1      | No stall |   0.009|
| YBL002W   | HTB2      | Stall    |   0.008|
| YER043C   | SAH1      | No stall |   0.008|
| YDR054C   | CDC34     | No stall |   0.008|
| YMR213W   | CEF1      | Stall    |   0.008|
| YOR245C   | DGA1      | No stall |   0.008|
| YBL068W   | PRS4      | No stall |   0.008|
| YBR021W   | FUR4      | No stall |   0.008|
| YDR498C   | SEC20     | No stall |   0.007|
| YDR150W   | NUM1      | No stall |   0.007|
| YDL042C   | SIR2      | No stall |   0.007|
| YBR228W   | SLX1      | No stall |   0.006|
| YOR141C   | ARP8      | Stall    |   0.006|
| YLR203C   | MSS51     | Stall    |   0.005|
| YDR085C   | AFR1      | Stall    |   0.005|
| YNL080C   | EOS1      | Stall    |   0.005|
| YPR163C   | TIF3      | No stall |   0.005|
| YEL038W   | UTR4      | No stall |   0.005|
| YBR102C   | EXO84     | No stall |   0.005|
| YIL005W   | EPS1      | Stall    |   0.004|
| YIL015W   | BAR1      | No stall |   0.004|
| YOL044W   | PEX15     | No stall |   0.004|
| YIL134W   | FLX1      | Stall    |   0.003|
| YBL045C   | COR1      | No stall |   0.002|
| YJL123C   | MTC1      | No stall |   0.002|
| YHR107C   | CDC12     | No stall |   0.001|
| YGL238W   | CSE1      | No stall |   0.001|
| YER015W   | FAA2      | No stall |   0.001|
| YHR092C   | HXT4      | No stall |   0.001|
| YDR055W   | PST1      | No stall |   0.000|
| YKR095W-A | PCC1      | No stall |  -0.001|
| YER118C   | SHO1      | No stall |  -0.001|
| YJR143C   | PMT4      | No stall |  -0.001|
| YCR018C   | SRD1      | Stall    |  -0.001|
| YDR177W   | UBC1      | No stall |  -0.001|
| YLR193C   | UPS1      | No stall |  -0.002|
| YGR091W   | PRP31     | Stall    |  -0.002|
| YDR283C   | GCN2      | Stall    |  -0.002|
| YCR092C   | MSH3      | Stall    |  -0.002|
| YAL034W-A | MTW1      | No stall |  -0.003|
| YJR051W   | OSM1      | No stall |  -0.003|
| YDR183W   | PLP1      | No stall |  -0.004|
| YBR283C   | SSH1      | No stall |  -0.004|
| YNL307C   | MCK1      | No stall |  -0.004|
| YLL001W   | DNM1      | No stall |  -0.004|
| YJL072C   | PSF2      | No stall |  -0.004|
| YDR169C   | STB3      | No stall |  -0.004|
| YBR203W   | COS111    | Stall    |  -0.005|
| YDR025W   | RPS11B    | Stall    |  -0.005|
| YNR009W   | NRM1      | No stall |  -0.006|
| YMR299C   | DYN3      | No stall |  -0.006|
| YOR334W   | MRS2      | Stall    |  -0.006|
| YOR005C   | DNL4      | Stall    |  -0.006|
| YGR103W   | NOP7      | Stall    |  -0.006|
| YER016W   | BIM1      | No stall |  -0.006|
| YJR053W   | BFA1      | No stall |  -0.007|
| YPL266W   | DIM1      | No stall |  -0.007|
| YGL133W   | ITC1      | Stall    |  -0.007|
| YOL137W   | BSC6      | No stall |  -0.007|
| YDR216W   | ADR1      | Stall    |  -0.007|
| YBR048W   | RPS11B    | Stall    |  -0.008|
| YLR238W   | FAR10     | Stall    |  -0.008|
| YDR528W   | HLR1      | Stall    |  -0.008|
| YKL219W   | COS9      | No stall |  -0.008|
| YBL032W   | HEK2      | No stall |  -0.008|
| YHR216W   | IMD2      | No stall |  -0.009|
| YIL159W   | BNR1      | Stall    |  -0.009|
| YOL040C   | RPS15     | No stall |  -0.009|
| YOL004W   | SIN3      | No stall |  -0.010|
| YBR166C   | TYR1      | No stall |  -0.010|
| YMR154C   | RIM13     | No stall |  -0.010|
| YAR002C-A | ERP1      | No stall |  -0.011|
| YMR066W   | SOV1      | Stall    |  -0.011|
| YFR028C   | CDC14     | Stall    |  -0.011|
| YGR271C-A | EFG1      | Stall    |  -0.011|
| YGR044C   | RME1      | No stall |  -0.011|
| YEL052W   | AFG1      | No stall |  -0.011|
| YMR100W   | MUB1      | Stall    |  -0.012|
| YGL066W   | SGF73     | Stall    |  -0.012|
| YBR046C   | ZTA1      | No stall |  -0.012|
| YDR532C   | KRE28     | No stall |  -0.013|
| YPR045C   | THP3      | Stall    |  -0.013|
| YDR448W   | ADA2      | No stall |  -0.013|
| YOL006C   | TOP1      | Stall    |  -0.013|
| YBL020W   | RFT1      | No stall |  -0.013|
| YMR270C   | RRN9      | Stall    |  -0.014|
| YJL192C   | SOP4      | No stall |  -0.014|
| YNL273W   | TOF1      | Stall    |  -0.014|
| YPR113W   | PIS1      | No stall |  -0.014|
| YBR094W   | PBY1      | No stall |  -0.015|
| YLR362W   | STE11     | No stall |  -0.015|
| YKR066C   | CCP1      | No stall |  -0.015|
| YKR096W   | ESL2      | Stall    |  -0.015|
| YLR348C   | DIC1      | No stall |  -0.015|
| YDL008W   | APC11     | No stall |  -0.015|
| YLR335W   | NUP2      | Stall    |  -0.015|
| YKL034W   | TUL1      | No stall |  -0.016|
| YPR104C   | FHL1      | Stall    |  -0.016|
| YNL189W   | SRP1      | No stall |  -0.016|
| YPR110C   | RPC40     | No stall |  -0.016|
| YNR016C   | ACC1      | No stall |  -0.016|
| YJL029C   | VPS53     | No stall |  -0.016|
| YER060W-A | FCY22     | No stall |  -0.016|
| YDR073W   | SNF11     | No stall |  -0.016|
| YOR113W   | AZF1      | No stall |  -0.017|
| YER032W   | FIR1      | Stall    |  -0.017|
| YIL018W   | RPL2A     | No stall |  -0.017|
| YKR022C   | NTR2      | No stall |  -0.017|
| YDR313C   | PIB1      | No stall |  -0.017|
| YBR058C   | UBP14     | No stall |  -0.017|
| YJR117W   | STE24     | No stall |  -0.018|
| YGR222W   | PET54     | No stall |  -0.018|
| YKL033W   | TTI1      | No stall |  -0.018|
| YDR259C   | YAP6      | No stall |  -0.018|
| YOR304W   | ISW2      | No stall |  -0.019|
| YJL177W   | RPL17B    | No stall |  -0.019|
| YGR178C   | PBP1      | No stall |  -0.019|
| YNR032W   | PPG1      | No stall |  -0.019|
| YMR088C   | VBA1      | No stall |  -0.019|
| YBR069C   | TAT1      | No stall |  -0.019|
| YGR138C   | TPO2      | Stall    |  -0.020|
| YGR084C   | MRP13     | No stall |  -0.021|
| YCR009C   | RVS161    | No stall |  -0.021|
| YDL082W   | RPL13A    | Stall    |  -0.021|
| YDR508C   | GNP1      | No stall |  -0.021|
| YKL120W   | OAC1      | No stall |  -0.021|
| YLL036C   | PRP19     | No stall |  -0.022|
| YGL207W   | SPT16     | No stall |  -0.022|
| YNL118C   | DCP2      | No stall |  -0.022|
| YDR251W   | PAM1      | Stall    |  -0.022|
| YDR137W   | RGP1      | No stall |  -0.022|
| YGR136W   | LSB1      | No stall |  -0.023|
| YGR252W   | GCN5      | No stall |  -0.023|
| YNL156C   | NSG2      | No stall |  -0.023|
| YHL039W   | EFM1      | No stall |  -0.024|
| YDR392W   | SPT3      | No stall |  -0.024|
| YHR007C   | ERG11     | No stall |  -0.025|
| YIL038C   | NOT3      | No stall |  -0.025|
| YFL023W   | BUD27     | Stall    |  -0.026|
| YGR189C   | CRH1      | No stall |  -0.026|
| YJL039C   | NUP192    | No stall |  -0.027|
| YLR078C   | BOS1      | No stall |  -0.027|
| YDL001W   | RMD1      | No stall |  -0.027|
| YDR057W   | YOS9      | No stall |  -0.027|
| YDR211W   | GCD6      | Stall    |  -0.028|
| YER006W   | NUG1      | Stall    |  -0.028|
| YOR073W   | SGO1      | Stall    |  -0.028|
| YNL250W   | RAD50     | No stall |  -0.028|
| YDR240C   | SNU56     | No stall |  -0.028|
| YGL022W   | STT3      | No stall |  -0.029|
| YBR227C   | MCX1      | No stall |  -0.029|
| YKR078W   | NA        | Stall    |  -0.030|
| YPL209C   | IPL1      | No stall |  -0.030|
| YIL138C   | TPM2      | No stall |  -0.030|
| YDL131W   | LYS21     | No stall |  -0.030|
| YDL040C   | NAT1      | No stall |  -0.031|
| YJL092W   | SRS2      | Stall    |  -0.031|
| YLR221C   | RSA3      | Stall    |  -0.031|
| YDL161W   | ENT1      | Stall    |  -0.031|
| YDR380W   | ARO10     | No stall |  -0.031|
| YGR040W   | KSS1      | No stall |  -0.032|
| YMR048W   | CSM3      | No stall |  -0.032|
| YER149C   | PEA2      | No stall |  -0.032|
| YNL247W   | NA        | No stall |  -0.032|
| YMR205C   | PFK2      | No stall |  -0.032|
| YPL122C   | TFB2      | Stall    |  -0.033|
| YLR023C   | IZH3      | No stall |  -0.033|
| YPR175W   | DPB2      | No stall |  -0.033|
| YKL172W   | EBP2      | Stall    |  -0.033|
| YOR290C   | SNF2      | Stall    |  -0.033|
| YIL035C   | CKA1      | Stall    |  -0.033|
| YKL024C   | URA6      | No stall |  -0.034|
| YKL145W   | RPT1      | No stall |  -0.034|
| YJR073C   | OPI3      | No stall |  -0.034|
| YAL026C   | DRS2      | No stall |  -0.034|
| YPR043W   | RPL43B    | No stall |  -0.035|
| YOL041C   | NOP12     | Stall    |  -0.035|
| YDR257C   | RKM4      | No stall |  -0.035|
| YML105C   | SEC65     | No stall |  -0.035|
| YNL280C   | ERG24     | No stall |  -0.035|
| YDR373W   | FRQ1      | No stall |  -0.035|
| YDR303C   | RSC3      | No stall |  -0.035|
| YKL060C   | FBA1      | No stall |  -0.035|
| YBR083W   | TEC1      | No stall |  -0.035|
| YNL032W   | SIW14     | No stall |  -0.036|
| YGR129W   | SYF2      | Stall    |  -0.036|
| YML118W   | NGL3      | No stall |  -0.036|
| YIL021W   | RPB3      | No stall |  -0.037|
| YLR276C   | DBP9      | Stall    |  -0.037|
| YDR130C   | FIN1      | No stall |  -0.037|
| YOR372C   | NDD1      | Stall    |  -0.038|
| YGR202C   | PCT1      | Stall    |  -0.038|
| YLR168C   | UPS2      | No stall |  -0.038|
| YER102W   | RPS8A     | Stall    |  -0.038|
| YLR165C   | PUS5      | No stall |  -0.038|
| YOL076W   | MDM20     | No stall |  -0.038|
| YNL291C   | MID1      | No stall |  -0.038|
| YGL196W   | DSD1      | No stall |  -0.038|
| YGR007W   | ECT1      | No stall |  -0.038|
| YBR192W   | RIM2      | No stall |  -0.039|
| YAR007C   | RFA1      | No stall |  -0.039|
| YDL160C   | DHH1      | No stall |  -0.039|
| YDL239C   | ADY3      | Stall    |  -0.040|
| YBL034C   | STU1      | No stall |  -0.040|
| YGL221C   | NIF3      | No stall |  -0.040|
| YOL058W   | ARG1      | No stall |  -0.040|
| YPL242C   | IQG1      | No stall |  -0.041|
| YLR154C   | RNH203    | No stall |  -0.042|
| YLR223C   | IFH1      | Stall    |  -0.042|
| YGR223C   | HSV2      | No stall |  -0.042|
| YAL042W   | ERV46     | No stall |  -0.043|
| YPL250C   | ATG41     | No stall |  -0.043|
| YNL246W   | VPS75     | No stall |  -0.043|
| YNL312W   | RFA2      | No stall |  -0.044|
| YJL131C   | AIM23     | Stall    |  -0.044|
| YJL156C   | SSY5      | No stall |  -0.044|
| YPR091C   | NVJ2      | Stall    |  -0.044|
| YPR019W   | MCM4      | No stall |  -0.044|
| YNR067C   | DSE4      | No stall |  -0.045|
| YJR066W   | TOR1      | No stall |  -0.045|
| YGR033C   | TIM21     | No stall |  -0.045|
| YER120W   | SCS2      | No stall |  -0.045|
| YGL194C   | HOS2      | No stall |  -0.046|
| YGL002W   | ERP6      | No stall |  -0.046|
| YFR047C   | BNA6      | No stall |  -0.046|
| YJR089W   | BIR1      | Stall    |  -0.046|
| YDR292C   | SRP101    | Stall    |  -0.046|
| YDL225W   | SHS1      | Stall    |  -0.046|
| YPL083C   | SEN54     | No stall |  -0.046|
| YIL016W   | SNL1      | Stall    |  -0.046|
| YGR102C   | GTF1      | No stall |  -0.047|
| YJL104W   | PAM16     | No stall |  -0.047|
| YGL065C   | ALG2      | No stall |  -0.047|
| YHR051W   | COX6      | No stall |  -0.047|
| YDR369C   | XRS2      | Stall    |  -0.047|
| YMR190C   | SGS1      | Stall    |  -0.047|
| YMR250W   | GAD1      | No stall |  -0.047|
| YLL061W   | MMP1      | No stall |  -0.048|
| YMR137C   | PSO2      | Stall    |  -0.048|
| YIR002C   | MPH1      | Stall    |  -0.048|
| YGR286C   | BIO2      | No stall |  -0.048|
| YPL084W   | BRO1      | No stall |  -0.048|
| YOL094C   | RFC4      | No stall |  -0.049|
| YPR023C   | EAF3      | No stall |  -0.049|
| YLR115W   | CFT2      | Stall    |  -0.049|
| YPL283C   | YRF1-7    | Stall    |  -0.050|
| YJR138W   | IML1      | Stall    |  -0.050|
| YJR014W   | TMA22     | Stall    |  -0.051|
| YLR128W   | DCN1      | No stall |  -0.052|
| YPL267W   | ACM1      | No stall |  -0.053|
| YJR059W   | PTK2      | No stall |  -0.053|
| YIL150C   | MCM10     | Stall    |  -0.053|
| YLR100W   | ERG27     | No stall |  -0.053|
| YDR062W   | LCB2      | No stall |  -0.053|
| YIL065C   | FIS1      | No stall |  -0.053|
| YDR524C   | AGE1      | No stall |  -0.054|
| YLL024C   | SSA2      | No stall |  -0.054|
| YOR033C   | EXO1      | Stall    |  -0.054|
| YJL088W   | ARG3      | No stall |  -0.054|
| YJL197W   | UBP12     | Stall    |  -0.055|
| YFR041C   | ERJ5      | No stall |  -0.055|
| YBR279W   | PAF1      | Stall    |  -0.055|
| YOR157C   | PUP1      | No stall |  -0.055|
| YMR283C   | RIT1      | No stall |  -0.055|
| YNL243W   | SLA2      | No stall |  -0.055|
| YKL213C   | DOA1      | No stall |  -0.055|
| YGL107C   | RMD9      | No stall |  -0.056|
| YDR322W   | MRPL35    | No stall |  -0.056|
| YGL164C   | YRB30     | No stall |  -0.056|
| YIL076W   | SEC28     | No stall |  -0.056|
| YBR059C   | AKL1      | No stall |  -0.056|
| YLR313C   | SPH1      | No stall |  -0.057|
| YJL035C   | TAD2      | No stall |  -0.057|
| YDL097C   | RPN6      | No stall |  -0.057|
| YER177W   | BMH1      | No stall |  -0.057|
| YOR118W   | RTC5      | No stall |  -0.057|
| YPR133C   | SPN1      | Stall    |  -0.057|
| YIL104C   | SHQ1      | No stall |  -0.057|
| YOR174W   | MED4      | No stall |  -0.058|
| YNL004W   | HRB1      | No stall |  -0.058|
| YDR231C   | COX20     | No stall |  -0.059|
| YFL037W   | TUB2      | No stall |  -0.059|
| YDL035C   | GPR1      | Stall    |  -0.060|
| YGL097W   | SRM1      | No stall |  -0.060|
| YCR037C   | PHO87     | Stall    |  -0.060|
| YKL042W   | SPC42     | Stall    |  -0.060|
| YKL189W   | HYM1      | No stall |  -0.060|
| YKL181W   | PRS1      | No stall |  -0.061|
| YDR376W   | ARH1      | No stall |  -0.061|
| YHR178W   | STB5      | Stall    |  -0.061|
| YJR043C   | POL32     | No stall |  -0.061|
| YLR195C   | NMT1      | No stall |  -0.061|
| YBR151W   | APD1      | No stall |  -0.061|
| YPL143W   | RPL33A    | No stall |  -0.062|
| YDR275W   | BSC2      | No stall |  -0.062|
| YOR323C   | PRO2      | No stall |  -0.062|
| YPL181W   | CTI6      | Stall    |  -0.062|
| YLR409C   | UTP21     | Stall    |  -0.062|
| YGL160W   | AIM14     | No stall |  -0.063|
| YOL145C   | CTR9      | Stall    |  -0.063|
| YBR070C   | ALG14     | No stall |  -0.063|
| YMR177W   | MMT1      | No stall |  -0.064|
| YMR054W   | STV1      | No stall |  -0.064|
| YGL123W   | RPS2      | No stall |  -0.064|
| YDL190C   | UFD2      | No stall |  -0.064|
| YDR108W   | TRS85     | No stall |  -0.064|
| YBR023C   | CHS3      | No stall |  -0.065|
| YGR285C   | ZUO1      | Stall    |  -0.065|
| YMR184W   | ADD37     | Stall    |  -0.065|
| YER072W   | VTC1      | No stall |  -0.065|
| YDL036C   | PUS9      | Stall    |  -0.065|
| YOL115W   | PAP2      | Stall    |  -0.066|
| YLR088W   | GAA1      | No stall |  -0.066|
| YMR075W   | RCO1      | Stall    |  -0.066|
| YNL131W   | TOM22     | No stall |  -0.066|
| YDR335W   | MSN5      | No stall |  -0.067|
| YLR201C   | COQ9      | No stall |  -0.067|
| YMR014W   | BUD22     | Stall    |  -0.067|
| YJL159W   | HSP150    | No stall |  -0.067|
| YLR007W   | NSE1      | No stall |  -0.067|
| YLR395C   | COX8      | No stall |  -0.068|
| YIL002C   | INP51     | No stall |  -0.069|
| YPR068C   | HOS1      | No stall |  -0.069|
| YPL045W   | VPS16     | No stall |  -0.070|
| YIL014W   | MNT3      | No stall |  -0.070|
| YOR074C   | CDC21     | No stall |  -0.071|
| YFR031C-A | RPL2A     | No stall |  -0.071|
| YOL039W   | RPP2A     | No stall |  -0.071|
| YDR213W   | UPC2      | Stall    |  -0.071|
| YML065W   | ORC1      | Stall    |  -0.072|
| YJL200C   | ACO2      | No stall |  -0.073|
| YOL005C   | RPB11     | No stall |  -0.073|
| YMR229C   | RRP5      | Stall    |  -0.073|
| YGL045W   | RIM8      | No stall |  -0.073|
| YNL157W   | IGO1      | No stall |  -0.073|
| YMR219W   | ESC1      | Stall    |  -0.073|
| YNL027W   | CRZ1      | Stall    |  -0.074|
| YGL067W   | NPY1      | No stall |  -0.074|
| YOL139C   | CDC33     | No stall |  -0.075|
| YCL016C   | DCC1      | Stall    |  -0.075|
| YDR163W   | CWC15     | No stall |  -0.075|
| YPR033C   | HTS1      | No stall |  -0.075|
| YDL055C   | PSA1      | No stall |  -0.076|
| YEL026W   | SNU13     | No stall |  -0.076|
| YDL174C   | DLD1      | Stall    |  -0.077|
| YGL092W   | NUP145    | No stall |  -0.077|
| YBR217W   | ATG12     | No stall |  -0.077|
| YBR245C   | ISW1      | No stall |  -0.078|
| YML071C   | COG8      | No stall |  -0.078|
| YDR295C   | HDA2      | Stall    |  -0.079|
| YGL061C   | DUO1      | No stall |  -0.079|
| YBR260C   | RGD1      | No stall |  -0.079|
| YPL154C   | PEP4      | No stall |  -0.079|
| YML100W   | TSL1      | No stall |  -0.080|
| YLR226W   | BUR2      | No stall |  -0.080|
| YGL195W   | GCN1      | No stall |  -0.080|
| YDL002C   | NHP10     | Stall    |  -0.081|
| YDR458C   | HEH2      | Stall    |  -0.081|
| YBR229C   | ROT2      | No stall |  -0.081|
| YHL048W   | COS8      | No stall |  -0.081|
| YLR389C   | STE23     | No stall |  -0.082|
| YDR093W   | DNF2      | Stall    |  -0.082|
| YPR061C   | JID1      | Stall    |  -0.082|
| YMR288W   | HSH155    | Stall    |  -0.083|
| YLR138W   | NHA1      | Stall    |  -0.083|
| YJL063C   | MRPL8     | No stall |  -0.084|
| YBR081C   | SPT7      | Stall    |  -0.084|
| YMR247C   | RKR1      | No stall |  -0.084|
| YLR266C   | PDR8      | Stall    |  -0.084|
| YMR108W   | ILV2      | No stall |  -0.085|
| YJL204C   | RCY1      | No stall |  -0.085|
| YLL043W   | FPS1      | Stall    |  -0.085|
| YNL059C   | ARP5      | Stall    |  -0.085|
| YER146W   | LSM5      | No stall |  -0.085|
| YDL135C   | RDI1      | No stall |  -0.085|
| YML070W   | DAK1      | No stall |  -0.085|
| YER031C   | YPT31     | No stall |  -0.086|
| YDL170W   | UGA3      | No stall |  -0.086|
| YGL068W   | MNP1      | No stall |  -0.086|
| YOL057W   | NA        | No stall |  -0.086|
| YGL143C   | MRF1      | Stall    |  -0.086|
| YDR159W   | SAC3      | No stall |  -0.087|
| YDL166C   | FAP7      | No stall |  -0.087|
| YDR298C   | ATP5      | No stall |  -0.087|
| YLR071C   | RGR1      | No stall |  -0.087|
| YDR288W   | NSE3      | No stall |  -0.088|
| YLL014W   | EMC6      | No stall |  -0.088|
| YER059W   | PCL6      | No stall |  -0.089|
| YNL208W   | NA        | No stall |  -0.089|
| YOR230W   | WTM1      | No stall |  -0.089|
| YDR190C   | RVB1      | No stall |  -0.089|
| YDR440W   | DOT1      | Stall    |  -0.089|
| YNL161W   | CBK1      | No stall |  -0.089|
| YDR429C   | TIF35     | No stall |  -0.090|
| YPR055W   | SEC8      | No stall |  -0.090|
| YOR151C   | RPB2      | No stall |  -0.091|
| YGL224C   | SDT1      | No stall |  -0.091|
| YDL116W   | NUP84     | No stall |  -0.091|
| YHR119W   | SET1      | Stall    |  -0.091|
| YJL168C   | SET2      | Stall    |  -0.091|
| YIL063C   | YRB2      | No stall |  -0.091|
| YHR190W   | ERG9      | No stall |  -0.092|
| YJR069C   | HAM1      | No stall |  -0.093|
| YOR286W   | RDL2      | No stall |  -0.093|
| YLR086W   | SMC4      | Stall    |  -0.093|
| YBR123C   | TFC1      | No stall |  -0.093|
| YBR154C   | RPB5      | No stall |  -0.094|
| YGR152C   | RSR1      | No stall |  -0.094|
| YDR099W   | BMH2      | No stall |  -0.094|
| YJR025C   | BNA1      | No stall |  -0.094|
| YER171W   | RAD3      | No stall |  -0.094|
| YMR022W   | UBC7      | No stall |  -0.095|
| YIR024C   | INA22     | No stall |  -0.095|
| YAL035W   | FUN12     | Stall    |  -0.095|
| YBR017C   | KAP104    | No stall |  -0.096|
| YAL031C   | GIP4      | Stall    |  -0.096|
| YDR364C   | CDC40     | Stall    |  -0.097|
| YEL058W   | PCM1      | No stall |  -0.097|
| YMR091C   | NPL6      | Stall    |  -0.097|
| YDR017C   | KCS1      | No stall |  -0.097|
| YHR076W   | PTC7      | No stall |  -0.098|
| YML062C   | MFT1      | Stall    |  -0.098|
| YNL267W   | PIK1      | No stall |  -0.099|
| YDR414C   | ERD1      | No stall |  -0.099|
| YKL130C   | SHE2      | No stall |  -0.099|
| YOR061W   | CKA2      | No stall |  -0.100|
| YKL215C   | OXP1      | No stall |  -0.100|
| YDL237W   | AIM6      | No stall |  -0.100|
| YIR003W   | AIM21     | Stall    |  -0.100|
| YBL009W   | ALK2      | No stall |  -0.100|
| YDL006W   | PTC1      | No stall |  -0.101|
| YPL138C   | SPP1      | Stall    |  -0.101|
| YKL122C   | SRP21     | Stall    |  -0.101|
| YDL113C   | ATG20     | No stall |  -0.101|
| YGR267C   | FOL2      | No stall |  -0.102|
| YLR250W   | SSP120    | No stall |  -0.102|
| YPL095C   | EEB1      | No stall |  -0.103|
| YDR113C   | PDS1      | No stall |  -0.103|
| YMR261C   | TPS3      | Stall    |  -0.103|
| YPL002C   | SNF8      | No stall |  -0.103|
| YHR151C   | MTC6      | No stall |  -0.103|
| YGL124C   | MON1      | No stall |  -0.103|
| YPR131C   | NAT3      | No stall |  -0.103|
| YGR240C   | PFK1      | No stall |  -0.103|
| YLR188W   | MDL1      | No stall |  -0.103|
| YHR005C   | GPA1      | No stall |  -0.103|
| YDR301W   | CFT1      | Stall    |  -0.104|
| YFR025C   | HIS2      | No stall |  -0.104|
| YGR167W   | CLC1      | No stall |  -0.104|
| YAL060W   | BDH1      | No stall |  -0.104|
| YML068W   | ITT1      | No stall |  -0.104|
| YBR068C   | BAP2      | No stall |  -0.104|
| YER164W   | CHD1      | Stall    |  -0.105|
| YGR054W   | NA        | Stall    |  -0.105|
| YGR163W   | GTR2      | No stall |  -0.105|
| YIL053W   | GPP1      | No stall |  -0.106|
| YGL125W   | MET13     | No stall |  -0.107|
| YDL095W   | PMT1      | No stall |  -0.107|
| YBR080C   | SEC18     | No stall |  -0.107|
| YOL038W   | PRE6      | No stall |  -0.108|
| YOR198C   | BFR1      | Stall    |  -0.108|
| YDR265W   | PEX10     | No stall |  -0.108|
| YDR058C   | TGL2      | No stall |  -0.108|
| YDR014W   | RAD61     | Stall    |  -0.109|
| YHR191C   | CTF8      | No stall |  -0.109|
| YCL032W   | STE50     | No stall |  -0.109|
| YHR063C   | PAN5      | No stall |  -0.109|
| YLR109W   | AHP1      | No stall |  -0.110|
| YJR006W   | POL31     | No stall |  -0.110|
| YDR482C   | CWC21     | No stall |  -0.110|
| YPL079W   | RPL21B    | No stall |  -0.110|
| YGR095C   | RRP46     | No stall |  -0.110|
| YPR129W   | SCD6      | No stall |  -0.111|
| YIL040W   | APQ12     | No stall |  -0.111|
| YJR131W   | MNS1      | No stall |  -0.111|
| YGL246C   | RAI1      | No stall |  -0.112|
| YLR180W   | SAM1      | No stall |  -0.112|
| YOR070C   | GYP1      | No stall |  -0.112|
| YER111C   | SWI4      | Stall    |  -0.112|
| YLR087C   | CSF1      | Stall    |  -0.112|
| YER148W   | SPT15     | No stall |  -0.112|
| YFL029C   | CAK1      | Stall    |  -0.112|
| YGR024C   | THG1      | No stall |  -0.113|
| YDR170C   | SEC7      | No stall |  -0.113|
| YPL058C   | PDR12     | No stall |  -0.113|
| YKL185W   | ASH1      | Stall    |  -0.113|
| YKL125W   | RRN3      | No stall |  -0.113|
| YPL176C   | TRE1      | Stall    |  -0.113|
| YCR033W   | SNT1      | Stall    |  -0.113|
| YPL167C   | REV3      | Stall    |  -0.113|
| YLR310C   | CDC25     | Stall    |  -0.113|
| YPR179C   | HDA3      | Stall    |  -0.114|
| YBL080C   | PET112    | No stall |  -0.114|
| YDL106C   | PHO2      | Stall    |  -0.114|
| YJR121W   | ATP2      | No stall |  -0.115|
| YOR347C   | PYK2      | No stall |  -0.115|
| YDL043C   | PRP11     | No stall |  -0.115|
| YGR258C   | RAD2      | Stall    |  -0.115|
| YDL104C   | QRI7      | No stall |  -0.115|
| YML004C   | GLO1      | No stall |  -0.117|
| YNL021W   | HDA1      | Stall    |  -0.117|
| YNL244C   | SUI1      | No stall |  -0.117|
| YOR350C   | MNE1      | No stall |  -0.117|
| YCL054W   | SPB1      | Stall    |  -0.117|
| YJL012C   | VTC4      | No stall |  -0.118|
| YGL040C   | HEM2      | No stall |  -0.118|
| YCR032W   | BPH1      | No stall |  -0.119|
| YFL039C   | ACT1      | No stall |  -0.119|
| YNL299W   | TRF5      | Stall    |  -0.119|
| YDR454C   | GUK1      | No stall |  -0.120|
| YNL323W   | LEM3      | Stall    |  -0.120|
| YNR028W   | CPR8      | No stall |  -0.120|
| YKL028W   | TFA1      | No stall |  -0.121|
| YGR254W   | ENO1      | No stall |  -0.121|
| YNL101W   | AVT4      | No stall |  -0.121|
| YJL115W   | ASF1      | No stall |  -0.121|
| YKL022C   | CDC16     | No stall |  -0.121|
| YNL083W   | SAL1      | No stall |  -0.121|
| YGR046W   | TAM41     | No stall |  -0.121|
| YLR450W   | HMG2      | No stall |  -0.121|
| YKR052C   | MRS4      | No stall |  -0.122|
| YKL116C   | PRR1      | No stall |  -0.122|
| YGL212W   | VAM7      | No stall |  -0.122|
| YBL103C   | RTG3      | Stall    |  -0.123|
| YHL019C   | APM2      | Stall    |  -0.123|
| YJL002C   | OST1      | No stall |  -0.123|
| YPL155C   | KIP2      | No stall |  -0.123|
| YDR088C   | SLU7      | Stall    |  -0.124|
| YLR357W   | RSC2      | Stall    |  -0.124|
| YPR140W   | TAZ1      | No stall |  -0.124|
| YNL304W   | YPT11     | Stall    |  -0.124|
| YPL240C   | HSP82     | Stall    |  -0.125|
| YFR001W   | LOC1      | Stall    |  -0.125|
| YPR035W   | GLN1      | Stall    |  -0.126|
| YMR192W   | GYL1      | Stall    |  -0.126|
| YDL240W   | LRG1      | No stall |  -0.126|
| YNR050C   | LYS9      | No stall |  -0.126|
| YDR043C   | NRG1      | No stall |  -0.127|
| YIL094C   | LYS12     | No stall |  -0.127|
| YDR264C   | AKR1      | No stall |  -0.127|
| YDR481C   | PHO8      | Stall    |  -0.127|
| YGR274C   | TAF1      | Stall    |  -0.127|
| YPR178W   | PRP4      | No stall |  -0.127|
| YEL018W   | EAF5      | No stall |  -0.128|
| YKL173W   | SNU114    | Stall    |  -0.128|
| YEL060C   | PRB1      | Stall    |  -0.128|
| YPL224C   | MMT2      | No stall |  -0.128|
| YNL102W   | POL1      | Stall    |  -0.128|
| YPR018W   | RLF2      | Stall    |  -0.129|
| YPL140C   | MKK2      | Stall    |  -0.129|
| YNL001W   | DOM34     | No stall |  -0.129|
| YDR127W   | ARO1      | No stall |  -0.129|
| YML031W   | NDC1      | No stall |  -0.129|
| YJR134C   | SGM1      | Stall    |  -0.129|
| YMR026C   | PEX12     | No stall |  -0.130|
| YKR003W   | OSH6      | No stall |  -0.130|
| YJL165C   | HAL5      | No stall |  -0.130|
| YDR477W   | SNF1      | No stall |  -0.130|
| YDR505C   | PSP1      | No stall |  -0.130|
| YKR095W   | MLP1      | No stall |  -0.130|
| YHR150W   | PEX28     | No stall |  -0.130|
| YHR142W   | CHS7      | No stall |  -0.130|
| YPL225W   | NA        | No stall |  -0.132|
| YLR344W   | RPL26A    | No stall |  -0.132|
| YNL301C   | RPL18B    | No stall |  -0.133|
| YLR314C   | CDC3      | No stall |  -0.133|
| YKL207W   | EMC3      | No stall |  -0.134|
| YPL011C   | TAF3      | No stall |  -0.134|
| YOR168W   | GLN4      | Stall    |  -0.134|
| YGL026C   | TRP5      | No stall |  -0.134|
| YDR261C   | EXG2      | No stall |  -0.134|
| YDR097C   | MSH6      | No stall |  -0.135|
| YBR204C   | LDH1      | No stall |  -0.135|
| YBR207W   | FTH1      | No stall |  -0.135|
| YIR011C   | STS1      | No stall |  -0.135|
| YDR311W   | TFB1      | No stall |  -0.135|
| YML046W   | PRP39     | No stall |  -0.135|
| YPR097W   | NA        | Stall    |  -0.136|
| YGL050W   | TYW3      | Stall    |  -0.136|
| YBR278W   | DPB3      | No stall |  -0.136|
| YFR024C-A | LSB3      | Stall    |  -0.137|
| YHR120W   | MSH1      | No stall |  -0.137|
| YGR122W   | NA        | No stall |  -0.137|
| YIL126W   | STH1      | Stall    |  -0.137|
| YFL013C   | IES1      | Stall    |  -0.138|
| YKL203C   | TOR2      | No stall |  -0.139|
| YFR053C   | HXK1      | No stall |  -0.139|
| YMR086W   | SEG1      | Stall    |  -0.139|
| YGR263C   | SAY1      | No stall |  -0.139|
| YGL245W   | GUS1      | No stall |  -0.139|
| YIL071C   | PCI8      | No stall |  -0.140|
| YCL011C   | GBP2      | Stall    |  -0.140|
| YMR078C   | CTF18     | Stall    |  -0.140|
| YKL165C   | MCD4      | No stall |  -0.141|
| YHR019C   | DED81     | No stall |  -0.141|
| YNL212W   | VID27     | No stall |  -0.141|
| YLR291C   | GCD7      | No stall |  -0.141|
| YOL117W   | RRI2      | No stall |  -0.141|
| YJR033C   | RAV1      | No stall |  -0.141|
| YCR031C   | RPS14A    | No stall |  -0.142|
| YLL015W   | BPT1      | No stall |  -0.142|
| YPL255W   | BBP1      | No stall |  -0.142|
| YDR182W   | CDC1      | Stall    |  -0.142|
| YOR166C   | SWT1      | No stall |  -0.143|
| YIL044C   | AGE2      | No stall |  -0.143|
| YHR196W   | UTP9      | No stall |  -0.143|
| YCR027C   | RHB1      | No stall |  -0.143|
| YMR237W   | BCH1      | No stall |  -0.143|
| YBL018C   | POP8      | No stall |  -0.143|
| YPL022W   | RAD1      | Stall    |  -0.143|
| YER005W   | YND1      | Stall    |  -0.144|
| YEL055C   | POL5      | Stall    |  -0.144|
| YBR233W-A | DAD3      | No stall |  -0.144|
| YDR164C   | SEC1      | Stall    |  -0.144|
| YGL257C   | MNT2      | No stall |  -0.144|
| YHR128W   | FUR1      | No stall |  -0.145|
| YER024W   | YAT2      | Stall    |  -0.145|
| YNL218W   | MGS1      | No stall |  -0.145|
| YER155C   | BEM2      | Stall    |  -0.146|
| YKL182W   | FAS1      | No stall |  -0.146|
| YGR130C   | NA        | No stall |  -0.146|
| YDL122W   | UBP1      | Stall    |  -0.147|
| YKR008W   | RSC4      | No stall |  -0.147|
| YHR024C   | MAS2      | No stall |  -0.147|
| YHL034C   | SBP1      | No stall |  -0.148|
| YML091C   | RPM2      | No stall |  -0.149|
| YPL074W   | YTA6      | No stall |  -0.149|
| YOR014W   | RTS1      | Stall    |  -0.149|
| YBR038W   | CHS2      | No stall |  -0.150|
| YDR438W   | THI74     | No stall |  -0.150|
| YNL072W   | RNH201    | Stall    |  -0.150|
| YNR017W   | TIM23     | No stall |  -0.150|
| YDR428C   | BNA7      | No stall |  -0.150|
| YOL097C   | WRS1      | No stall |  -0.151|
| YKR062W   | TFA2      | Stall    |  -0.151|
| YBR143C   | SUP45     | No stall |  -0.151|
| YKR100C   | SKG1      | No stall |  -0.151|
| YNL070W   | TOM7      | No stall |  -0.151|
| YML115C   | VAN1      | No stall |  -0.152|
| YML085C   | TUB1      | No stall |  -0.152|
| YJL187C   | SWE1      | Stall    |  -0.153|
| YIL085C   | KTR7      | No stall |  -0.153|
| YLR207W   | HRD3      | No stall |  -0.153|
| YBR087W   | RFC5      | No stall |  -0.154|
| YJL213W   | NA        | No stall |  -0.154|
| YLR273C   | PIG1      | No stall |  -0.154|
| YIL041W   | GVP36     | No stall |  -0.154|
| YER020W   | GPA2      | No stall |  -0.155|
| YKL168C   | KKQ8      | No stall |  -0.155|
| YLR067C   | PET309    | No stall |  -0.155|
| YGL252C   | RTG2      | No stall |  -0.155|
| YOR261C   | RPN8      | No stall |  -0.156|
| YDL045C   | FAD1      | No stall |  -0.156|
| YGL203C   | KEX1      | Stall    |  -0.157|
| YDR394W   | RPT3      | No stall |  -0.157|
| YPR118W   | MRI1      | No stall |  -0.158|
| YGR123C   | PPT1      | No stall |  -0.158|
| YML074C   | FPR3      | Stall    |  -0.158|
| YLR153C   | ACS2      | No stall |  -0.158|
| YML088W   | UFO1      | Stall    |  -0.159|
| YDR372C   | VPS74     | No stall |  -0.159|
| YLR179C   | NA        | No stall |  -0.159|
| YPR021C   | AGC1      | Stall    |  -0.159|
| YLR321C   | SFH1      | Stall    |  -0.159|
| YDR235W   | PRP42     | Stall    |  -0.160|
| YOL002C   | IZH2      | Stall    |  -0.160|
| YKR019C   | IRS4      | Stall    |  -0.160|
| YNL294C   | RIM21     | Stall    |  -0.160|
| YGR104C   | SRB5      | No stall |  -0.160|
| YBL069W   | AST1      | No stall |  -0.161|
| YPR173C   | VPS4      | No stall |  -0.162|
| YGL100W   | SEH1      | No stall |  -0.162|
| YDL185W   | VMA1      | No stall |  -0.162|
| YOR324C   | FRT1      | No stall |  -0.162|
| YKL049C   | CSE4      | Stall    |  -0.162|
| YIL112W   | HOS4      | Stall    |  -0.163|
| YNL038W   | GPI15     | No stall |  -0.163|
| YOL013C   | HRD1      | Stall    |  -0.163|
| YLR015W   | BRE2      | Stall    |  -0.163|
| YDL182W   | LYS20     | No stall |  -0.163|
| YMR162C   | DNF3      | No stall |  -0.164|
| YDR462W   | MRPL28    | No stall |  -0.164|
| YDR122W   | KIN1      | Stall    |  -0.164|
| YHR031C   | RRM3      | No stall |  -0.165|
| YOL011W   | PLB3      | No stall |  -0.166|
| YBL008W   | HIR1      | No stall |  -0.166|
| YPR101W   | SNT309    | Stall    |  -0.166|
| YCL038C   | ATG22     | No stall |  -0.167|
| YIL074C   | SER33     | No stall |  -0.167|
| YJR140C   | HIR3      | No stall |  -0.168|
| YGR112W   | SHY1      | No stall |  -0.168|
| YMR225C   | MRPL44    | No stall |  -0.168|
| YPL254W   | HFI1      | No stall |  -0.168|
| YMR286W   | MRPL33    | No stall |  -0.168|
| YKL157W   | APE2      | No stall |  -0.168|
| YOR191W   | ULS1      | Stall    |  -0.168|
| YDR004W   | RAD57     | No stall |  -0.168|
| YKL170W   | MRPL38    | No stall |  -0.169|
| YGL112C   | TAF6      | No stall |  -0.169|
| YJR058C   | APS2      | No stall |  -0.169|
| YGR246C   | BRF1      | Stall    |  -0.169|
| YCR002C   | CDC10     | No stall |  -0.170|
| YDR300C   | PRO1      | No stall |  -0.171|
| YOL148C   | SPT20     | Stall    |  -0.171|
| YLR275W   | SMD2      | No stall |  -0.171|
| YDR245W   | MNN10     | No stall |  -0.171|
| YML117W   | NAB6      | Stall    |  -0.171|
| YMR268C   | PRP24     | No stall |  -0.172|
| YBR273C   | UBX7      | No stall |  -0.172|
| YCL057W   | PRD1      | Stall    |  -0.172|
| YBL017C   | PEP1      | Stall    |  -0.172|
| YOL120C   | RPL18B    | No stall |  -0.172|
| YGL001C   | ERG26     | No stall |  -0.172|
| YGR262C   | BUD32     | No stall |  -0.172|
| YML097C   | VPS9      | No stall |  -0.173|
| YJR013W   | GPI14     | No stall |  -0.173|
| YML069W   | POB3      | No stall |  -0.173|
| YGR036C   | CAX4      | No stall |  -0.173|
| YBL058W   | SHP1      | No stall |  -0.174|
| YDR141C   | DOP1      | Stall    |  -0.174|
| YJR102C   | VPS25     | No stall |  -0.174|
| YIL131C   | FKH1      | No stall |  -0.175|
| YDR142C   | PEX7      | No stall |  -0.175|
| YMR012W   | CLU1      | No stall |  -0.175|
| YKL096W-A | CWP2      | No stall |  -0.175|
| YGR196C   | FYV8      | No stall |  -0.175|
| YKL210W   | UBA1      | No stall |  -0.176|
| YBL072C   | RPS8A     | Stall    |  -0.176|
| YOR287C   | RRP36     | Stall    |  -0.177|
| YBR275C   | RIF1      | Stall    |  -0.177|
| YJL143W   | TIM17     | No stall |  -0.177|
| YGL016W   | KAP122    | No stall |  -0.177|
| YHR027C   | RPN1      | No stall |  -0.178|
| YBR216C   | YBP1      | No stall |  -0.178|
| YNL073W   | MSK1      | No stall |  -0.178|
| YDR363W-A | SEM1      | No stall |  -0.178|
| YEL013W   | VAC8      | No stall |  -0.178|
| YDL134C   | PPH21     | No stall |  -0.178|
| YGR218W   | CRM1      | No stall |  -0.178|
| YAL016W   | TPD3      | No stall |  -0.178|
| YBR145W   | ADH5      | No stall |  -0.178|
| YGR098C   | ESP1      | Stall    |  -0.179|
| YOL049W   | GSH2      | No stall |  -0.179|
| YIR025W   | MND2      | Stall    |  -0.179|
| YDR315C   | IPK1      | No stall |  -0.180|
| YBR290W   | BSD2      | No stall |  -0.180|
| YNL325C   | FIG4      | No stall |  -0.180|
| YMR230W   | RPS10B    | No stall |  -0.180|
| YLR029C   | RPL15A    | Stall    |  -0.180|
| YBR159W   | IFA38     | No stall |  -0.180|
| YIL043C   | CBR1      | No stall |  -0.181|
| YNR053C   | NOG2      | Stall    |  -0.181|
| YLL039C   | UBI4      | No stall |  -0.181|
| YLR277C   | YSH1      | No stall |  -0.181|
| YOL088C   | MPD2      | No stall |  -0.182|
| YGR034W   | RPL26B    | No stall |  -0.182|
| YDR217C   | RAD9      | Stall    |  -0.182|
| YNL029C   | KTR5      | Stall    |  -0.183|
| YOL129W   | VPS68     | No stall |  -0.183|
| YNL311C   | SKP2      | No stall |  -0.183|
| YPR081C   | GRS2      | No stall |  -0.183|
| YCR023C   | NA        | No stall |  -0.183|
| YFR043C   | IRC6      | No stall |  -0.184|
| YCR017C   | CWH43     | No stall |  -0.184|
| YIL157C   | COA1      | No stall |  -0.184|
| YOR327C   | SNC2      | No stall |  -0.184|
| YJR024C   | MDE1      | No stall |  -0.184|
| YDR304C   | CPR5      | No stall |  -0.185|
| YGL175C   | SAE2      | Stall    |  -0.186|
| YHR166C   | CDC23     | No stall |  -0.186|
| YPR176C   | BET2      | No stall |  -0.186|
| YMR125W   | STO1      | No stall |  -0.186|
| YMR235C   | RNA1      | No stall |  -0.186|
| YGL191W   | COX13     | No stall |  -0.186|
| YKL073W   | LHS1      | Stall    |  -0.187|
| YPL173W   | MRPL40    | No stall |  -0.187|
| YOR216C   | RUD3      | Stall    |  -0.187|
| YNL188W   | KAR1      | No stall |  -0.187|
| YGR186W   | TFG1      | Stall    |  -0.188|
| YLR247C   | IRC20     | Stall    |  -0.188|
| YDR331W   | GPI8      | No stall |  -0.188|
| YLR440C   | SEC39     | No stall |  -0.188|
| YBR179C   | FZO1      | Stall    |  -0.189|
| YPL029W   | SUV3      | No stall |  -0.189|
| YDL193W   | NUS1      | No stall |  -0.189|
| YPL170W   | DAP1      | No stall |  -0.189|
| YML029W   | USA1      | No stall |  -0.189|
| YOR172W   | YRM1      | Stall    |  -0.190|
| YOR089C   | VPS21     | No stall |  -0.190|
| YNL239W   | LAP3      | No stall |  -0.190|
| YHR129C   | ARP1      | No stall |  -0.190|
| YHR152W   | SPO12     | Stall    |  -0.191|
| YNL216W   | RAP1      | No stall |  -0.191|
| YHR106W   | TRR2      | No stall |  -0.191|
| YDR180W   | SCC2      | Stall    |  -0.192|
| YMR197C   | VTI1      | No stall |  -0.192|
| YOR374W   | ALD4      | No stall |  -0.192|
| YJR005W   | APL1      | No stall |  -0.193|
| YNL242W   | ATG2      | No stall |  -0.193|
| YGR169C   | PUS6      | No stall |  -0.194|
| YGR295C   | COS6      | No stall |  -0.194|
| YGR184C   | UBR1      | Stall    |  -0.194|
| YEL056W   | HAT2      | No stall |  -0.194|
| YOL054W   | PSH1      | No stall |  -0.195|
| YFL016C   | MDJ1      | No stall |  -0.195|
| YMR094W   | CTF13     | Stall    |  -0.196|
| YML126C   | ERG13     | No stall |  -0.196|
| YMR200W   | ROT1      | No stall |  -0.196|
| YOL103W   | ITR2      | No stall |  -0.197|
| YBL046W   | PSY4      | Stall    |  -0.197|
| YLL006W   | MMM1      | No stall |  -0.197|
| YDL117W   | CYK3      | Stall    |  -0.197|
| YBL078C   | ATG8      | No stall |  -0.197|
| YBR121C   | GRS1      | Stall    |  -0.198|
| YDR430C   | CYM1      | No stall |  -0.198|
| YMR161W   | HLJ1      | No stall |  -0.198|
| YFR032C-A | RPL29     | No stall |  -0.198|
| YAR018C   | KIN3      | No stall |  -0.198|
| YGR076C   | MRPL25    | No stall |  -0.198|
| YAR008W   | SEN34     | No stall |  -0.198|
| YBR251W   | MRPS5     | No stall |  -0.198|
| YLR032W   | RAD5      | No stall |  -0.198|
| YDL176W   | NA        | No stall |  -0.199|
| YGL035C   | MIG1      | Stall    |  -0.199|
| YLR435W   | TSR2      | Stall    |  -0.199|
| YFL047W   | RGD2      | No stall |  -0.199|
| YHR174W   | ENO2      | No stall |  -0.200|
| YKL166C   | TPK3      | No stall |  -0.200|
| YDR236C   | FMN1      | No stall |  -0.201|
| YNL127W   | FAR11     | No stall |  -0.201|
| YMR294W   | JNM1      | No stall |  -0.202|
| YBL061C   | SKT5      | Stall    |  -0.202|
| YPL072W   | UBP16     | Stall    |  -0.202|
| YFR030W   | MET10     | No stall |  -0.202|
| YDR107C   | TMN2      | No stall |  -0.203|
| YJL036W   | SNX4      | No stall |  -0.203|
| YKL079W   | SMY1      | No stall |  -0.204|
| YDL209C   | CWC2      | No stall |  -0.204|
| YLR438C-A | LSM3      | No stall |  -0.205|
| YKL048C   | ELM1      | Stall    |  -0.205|
| YBL006C   | LDB7      | Stall    |  -0.205|
| YMR211W   | DML1      | No stall |  -0.205|
| YDL136W   | RPL35A    | Stall    |  -0.205|
| YIR038C   | GTT1      | No stall |  -0.205|
| YIR005W   | IST3      | No stall |  -0.206|
| YOR090C   | PTC5      | No stall |  -0.206|
| YOR332W   | VMA4      | No stall |  -0.206|
| YGL020C   | GET1      | No stall |  -0.206|
| YGR013W   | SNU71     | No stall |  -0.207|
| YMR220W   | ERG8      | No stall |  -0.207|
| YHR050W   | SMF2      | No stall |  -0.208|
| YOR035C   | SHE4      | No stall |  -0.208|
| YBR196C   | PGI1      | No stall |  -0.208|
| YCL012C   | NA        | No stall |  -0.209|
| YBR236C   | ABD1      | No stall |  -0.209|
| YPL211W   | NIP7      | No stall |  -0.210|
| YMR295C   | NA        | Stall    |  -0.210|
| YOL078W   | AVO1      | No stall |  -0.210|
| YLR387C   | REH1      | Stall    |  -0.210|
| YHR110W   | ERP5      | No stall |  -0.210|
| YDL155W   | CLB3      | No stall |  -0.211|
| YJL171C   | TOH1      | No stall |  -0.211|
| YPR028W   | YOP1      | No stall |  -0.211|
| YLR433C   | CNA1      | No stall |  -0.211|
| YBR049C   | REB1      | Stall    |  -0.212|
| YPL024W   | RMI1      | No stall |  -0.212|
| YMR028W   | TAP42     | No stall |  -0.213|
| YPL268W   | PLC1      | Stall    |  -0.213|
| YBR002C   | RER2      | No stall |  -0.213|
| YNL272C   | SEC2      | Stall    |  -0.213|
| YOR095C   | RKI1      | No stall |  -0.214|
| YML076C   | WAR1      | Stall    |  -0.214|
| YDL132W   | CDC53     | No stall |  -0.214|
| YML038C   | YMD8      | No stall |  -0.216|
| YIL007C   | NAS2      | No stall |  -0.217|
| YIL022W   | TIM44     | Stall    |  -0.217|
| YPR167C   | MET16     | No stall |  -0.217|
| YER048C   | CAJ1      | No stall |  -0.217|
| YOR270C   | VPH1      | No stall |  -0.218|
| YAL024C   | LTE1      | No stall |  -0.218|
| YMR011W   | HXT2      | No stall |  -0.218|
| YLR319C   | BUD6      | No stall |  -0.219|
| YKL149C   | DBR1      | No stall |  -0.219|
| YPR187W   | RPO26     | No stall |  -0.219|
| YHR156C   | LIN1      | Stall    |  -0.219|
| YGR261C   | APL6      | Stall    |  -0.219|
| YML061C   | PIF1      | Stall    |  -0.219|
| YER053C   | PIC2      | No stall |  -0.219|
| YDR096W   | GIS1      | Stall    |  -0.219|
| YDR330W   | UBX5      | No stall |  -0.220|
| YGR209C   | TRX2      | No stall |  -0.220|
| YDR243C   | PRP28     | No stall |  -0.220|
| YKL143W   | LTV1      | Stall    |  -0.220|
| YIR010W   | DSN1      | Stall    |  -0.220|
| YDR080W   | VPS41     | Stall    |  -0.221|
| YNL098C   | RAS2      | No stall |  -0.222|
| YDL159W   | STE7      | Stall    |  -0.222|
| YLL004W   | ORC3      | No stall |  -0.223|
| YFR007W   | YFH7      | No stall |  -0.223|
| YMR264W   | CUE1      | No stall |  -0.223|
| YNR037C   | RSM19     | No stall |  -0.223|
| YLR289W   | GUF1      | Stall    |  -0.223|
| YER136W   | GDI1      | No stall |  -0.223|
| YKL039W   | PTM1      | No stall |  -0.223|
| YJR126C   | VPS70     | No stall |  -0.224|
| YOR234C   | RPL33B    | No stall |  -0.224|
| YNL290W   | RFC3      | No stall |  -0.224|
| YMR210W   | MGL2      | Stall    |  -0.225|
| YKR001C   | VPS1      | No stall |  -0.225|
| YPL001W   | HAT1      | No stall |  -0.225|
| YOL026C   | MIM1      | No stall |  -0.225|
| YMR020W   | FMS1      | No stall |  -0.225|
| YBL023C   | MCM2      | Stall    |  -0.226|
| YGL038C   | OCH1      | No stall |  -0.226|
| YBL041W   | PRE7      | No stall |  -0.226|
| YOR326W   | MYO2      | Stall    |  -0.226|
| YLR117C   | CLF1      | Stall    |  -0.226|
| YBR098W   | MMS4      | No stall |  -0.227|
| YFR014C   | CMK1      | No stall |  -0.227|
| YGR232W   | NAS6      | No stall |  -0.227|
| YCR030C   | SYP1      | Stall    |  -0.228|
| YML028W   | TSA1      | No stall |  -0.228|
| YPL203W   | TPK2      | No stall |  -0.228|
| YDR325W   | YCG1      | Stall    |  -0.228|
| YLR399C   | BDF1      | Stall    |  -0.229|
| YLL011W   | SOF1      | Stall    |  -0.229|
| YKL087C   | CYT2      | No stall |  -0.229|
| YGL142C   | GPI10     | No stall |  -0.229|
| YNL178W   | RPS3      | No stall |  -0.230|
| YML104C   | MDM1      | No stall |  -0.230|
| YJL194W   | CDC6      | Stall    |  -0.231|
| YJR031C   | GEA1      | No stall |  -0.231|
| YLR417W   | VPS36     | No stall |  -0.231|
| YCR094W   | CDC50     | No stall |  -0.231|
| YGR157W   | CHO2      | No stall |  -0.231|
| YJL030W   | MAD2      | No stall |  -0.231|
| YML094W   | GIM5      | No stall |  -0.231|
| YIL137C   | TMA108    | No stall |  -0.232|
| YLR392C   | ART10     | No stall |  -0.232|
| YPL215W   | CBP3      | No stall |  -0.233|
| YMR224C   | MRE11     | Stall    |  -0.233|
| YKR014C   | YPT52     | No stall |  -0.233|
| YER019C-A | SBH2      | No stall |  -0.234|
| YKL069W   | NA        | No stall |  -0.235|
| YPL019C   | VTC3      | Stall    |  -0.235|
| YLR028C   | ADE16     | No stall |  -0.235|
| YAL011W   | SWC3      | Stall    |  -0.235|
| YNL051W   | COG5      | No stall |  -0.235|
| YDL194W   | SNF3      | No stall |  -0.236|
| YBR257W   | POP4      | Stall    |  -0.236|
| YGR028W   | MSP1      | No stall |  -0.236|
| YKL154W   | SRP102    | No stall |  -0.236|
| YDR003W   | RCR2      | No stall |  -0.236|
| YOR264W   | DSE3      | No stall |  -0.236|
| YMR041C   | ARA2      | No stall |  -0.236|
| YJL008C   | CCT8      | No stall |  -0.237|
| YIL144W   | NDC80     | No stall |  -0.237|
| YLR039C   | RIC1      | No stall |  -0.237|
| YKL208W   | CBT1      | No stall |  -0.237|
| YJL190C   | RPS22A    | No stall |  -0.237|
| YOR026W   | BUB3      | No stall |  -0.238|
| YPL004C   | LSP1      | No stall |  -0.238|
| YDR044W   | HEM13     | No stall |  -0.238|
| YMR226C   | NA        | No stall |  -0.238|
| YOL082W   | ATG19     | No stall |  -0.239|
| YLR130C   | ZRT2      | No stall |  -0.239|
| YKR063C   | LAS1      | Stall    |  -0.239|
| YOR285W   | RDL1      | No stall |  -0.240|
| YBR140C   | IRA1      | No stall |  -0.240|
| YDR191W   | HST4      | Stall    |  -0.240|
| YDL219W   | DTD1      | No stall |  -0.240|
| YGL173C   | XRN1      | Stall    |  -0.240|
| YDR468C   | TLG1      | No stall |  -0.240|
| YKL011C   | CCE1      | No stall |  -0.240|
| YDR074W   | TPS2      | No stall |  -0.241|
| YGL077C   | HNM1      | No stall |  -0.241|
| YBR274W   | CHK1      | No stall |  -0.241|
| YOR171C   | LCB4      | No stall |  -0.242|
| YPL235W   | RVB2      | No stall |  -0.242|
| YDL110C   | TMA17     | No stall |  -0.243|
| YJL196C   | ELO1      | No stall |  -0.243|
| YJL046W   | AIM22     | No stall |  -0.243|
| YFL038C   | YPT1      | No stall |  -0.243|
| YJL054W   | TIM54     | No stall |  -0.244|
| YKL134C   | OCT1      | No stall |  -0.245|
| YML127W   | RSC9      | No stall |  -0.245|
| YLR429W   | CRN1      | No stall |  -0.245|
| YGR206W   | MVB12     | No stall |  -0.245|
| YGL228W   | SHE10     | Stall    |  -0.246|
| YEL002C   | WBP1      | No stall |  -0.246|
| YCL061C   | MRC1      | Stall    |  -0.246|
| YGR276C   | RNH70     | Stall    |  -0.246|
| YDR470C   | UGO1      | No stall |  -0.247|
| YER073W   | ALD5      | No stall |  -0.247|
| YLR182W   | SWI6      | No stall |  -0.248|
| YDR485C   | VPS72     | Stall    |  -0.248|
| YPR088C   | SRP54     | No stall |  -0.249|
| YGL095C   | VPS45     | No stall |  -0.249|
| YFL027C   | GYP8      | No stall |  -0.250|
| YJL186W   | MNN5      | No stall |  -0.250|
| YDL067C   | COX9      | No stall |  -0.250|
| YGR278W   | CWC22     | No stall |  -0.251|
| YHR036W   | BRL1      | No stall |  -0.251|
| YDR473C   | PRP3      | Stall    |  -0.251|
| YGL141W   | HUL5      | Stall    |  -0.251|
| YML106W   | URA5      | No stall |  -0.252|
| YDR503C   | LPP1      | Stall    |  -0.252|
| YMR157C   | AIM36     | No stall |  -0.252|
| YDL053C   | PBP4      | No stall |  -0.252|
| YDR212W   | TCP1      | Stall    |  -0.252|
| YGR100W   | MDR1      | No stall |  -0.252|
| YDL022W   | GPD1      | No stall |  -0.253|
| YDR490C   | PKH1      | No stall |  -0.253|
| YKL196C   | YKT6      | No stall |  -0.253|
| YGR171C   | MSM1      | No stall |  -0.253|
| YLR386W   | VAC14     | No stall |  -0.253|
| YKL062W   | MSN4      | No stall |  -0.253|
| YLR099C   | ICT1      | No stall |  -0.254|
| YDL212W   | SHR3      | No stall |  -0.254|
| YIL008W   | URM1      | No stall |  -0.255|
| YDL019C   | OSH2      | Stall    |  -0.255|
| YBR157C   | ICS2      | No stall |  -0.255|
| YLR216C   | CPR6      | No stall |  -0.255|
| YDR471W   | RPL27B    | Stall    |  -0.255|
| YBR200W   | BEM1      | No stall |  -0.256|
| YDR478W   | SNM1      | Stall    |  -0.256|
| YMR183C   | SSO2      | Stall    |  -0.256|
| YAL001C   | TFC3      | Stall    |  -0.256|
| YOL056W   | GPM3      | No stall |  -0.257|
| YDL173W   | PAR32     | Stall    |  -0.258|
| YLR421C   | RPN13     | No stall |  -0.258|
| YBR254C   | TRS20     | No stall |  -0.258|
| YKL112W   | ABF1      | Stall    |  -0.258|
| YER183C   | FAU1      | No stall |  -0.259|
| YKL035W   | UGP1      | No stall |  -0.259|
| YOR281C   | PLP2      | No stall |  -0.259|
| YJR065C   | ARP3      | No stall |  -0.259|
| YDR435C   | PPM1      | No stall |  -0.260|
| YBR042C   | CST26     | No stall |  -0.260|
| YOL018C   | TLG2      | No stall |  -0.261|
| YLR119W   | SRN2      | Stall    |  -0.261|
| YDR013W   | PSF1      | No stall |  -0.261|
| YDL102W   | POL3      | Stall    |  -0.262|
| YJL137C   | GLG2      | No stall |  -0.263|
| YKL089W   | MIF2      | Stall    |  -0.263|
| YMR289W   | ABZ2      | No stall |  -0.263|
| YGL155W   | CDC43     | No stall |  -0.264|
| YDR318W   | MCM21     | Stall    |  -0.264|
| YOR371C   | GPB1      | No stall |  -0.264|
| YCL031C   | RRP7      | Stall    |  -0.264|
| YIL036W   | CST6      | No stall |  -0.264|
| YEL024W   | RIP1      | No stall |  -0.265|
| YNR020C   | ATP23     | No stall |  -0.265|
| YGR140W   | CBF2      | Stall    |  -0.265|
| YPL015C   | HST2      | No stall |  -0.265|
| YDL234C   | GYP7      | No stall |  -0.265|
| YJL034W   | KAR2      | No stall |  -0.265|
| YPR105C   | COG4      | No stall |  -0.265|
| YDL101C   | DUN1      | No stall |  -0.266|
| YNR031C   | SSK2      | No stall |  -0.266|
| YJR125C   | ENT3      | Stall    |  -0.266|
| YPR049C   | ATG11     | No stall |  -0.267|
| YPL082C   | MOT1      | Stall    |  -0.267|
| YOR346W   | REV1      | No stall |  -0.267|
| YML016C   | PPZ1      | No stall |  -0.267|
| YGR217W   | CCH1      | No stall |  -0.267|
| YGL134W   | PCL10     | No stall |  -0.267|
| YGR005C   | TFG2      | Stall    |  -0.267|
| YMR105C   | PGM2      | No stall |  -0.268|
| YMR109W   | MYO5      | Stall    |  -0.268|
| YOR356W   | CIR2      | No stall |  -0.268|
| YJL176C   | SWI3      | No stall |  -0.268|
| YIL069C   | RPS24A    | Stall    |  -0.269|
| YDR105C   | TMS1      | No stall |  -0.269|
| YMR207C   | HFA1      | No stall |  -0.269|
| YDR334W   | SWR1      | Stall    |  -0.269|
| YDR457W   | TOM1      | Stall    |  -0.269|
| YHR113W   | APE4      | No stall |  -0.269|
| YCL017C   | NFS1      | No stall |  -0.269|
| YLR298C   | YHC1      | Stall    |  -0.270|
| YHR034C   | PIH1      | No stall |  -0.270|
| YOR274W   | MOD5      | No stall |  -0.270|
| YPR108W   | RPN7      | No stall |  -0.270|
| YDL093W   | PMT5      | No stall |  -0.270|
| YPL196W   | OXR1      | No stall |  -0.271|
| YBL067C   | UBP13     | Stall    |  -0.271|
| YMR240C   | CUS1      | Stall    |  -0.271|
| YNL227C   | JJJ1      | Stall    |  -0.271|
| YDR379C-A | SDH6      | No stall |  -0.271|
| YER030W   | CHZ1      | Stall    |  -0.271|
| YDR002W   | YRB1      | No stall |  -0.272|
| YLR021W   | IRC25     | No stall |  -0.272|
| YCR011C   | ADP1      | No stall |  -0.272|
| YLR163C   | MAS1      | No stall |  -0.272|
| YLR324W   | PEX30     | Stall    |  -0.272|
| YDL125C   | HNT1      | No stall |  -0.273|
| YPL237W   | SUI3      | Stall    |  -0.273|
| YKR006C   | MRPL13    | Stall    |  -0.274|
| YBR244W   | GPX2      | No stall |  -0.274|
| YDL145C   | COP1      | No stall |  -0.274|
| YPR138C   | MEP3      | No stall |  -0.274|
| YBR095C   | RXT2      | No stall |  -0.274|
| YGL179C   | TOS3      | No stall |  -0.275|
| YML072C   | TCB3      | Stall    |  -0.275|
| YPR152C   | URN1      | No stall |  -0.275|
| YGL084C   | GUP1      | No stall |  -0.275|
| YOR069W   | VPS5      | No stall |  -0.275|
| YMR003W   | AIM34     | No stall |  -0.275|
| YHR026W   | VMA16     | No stall |  -0.276|
| YER166W   | DNF1      | Stall    |  -0.276|
| YDR427W   | RPN9      | No stall |  -0.277|
| YDL078C   | MDH3      | No stall |  -0.277|
| YNR033W   | ABZ1      | No stall |  -0.277|
| YOL101C   | IZH4      | No stall |  -0.277|
| YFR009W   | GCN20     | No stall |  -0.278|
| YNL100W   | MIC27     | No stall |  -0.278|
| YPL156C   | PRM4      | No stall |  -0.278|
| YJL183W   | MNN11     | No stall |  -0.278|
| YDR475C   | JIP4      | Stall    |  -0.278|
| YGL137W   | SEC27     | No stall |  -0.278|
| YLR118C   | NA        | No stall |  -0.279|
| YLR097C   | HRT3      | No stall |  -0.279|
| YKL101W   | HSL1      | No stall |  -0.279|
| YER026C   | CHO1      | No stall |  -0.279|
| YGR187C   | HGH1      | No stall |  -0.279|
| YOL033W   | MSE1      | No stall |  -0.280|
| YLR082C   | SRL2      | No stall |  -0.280|
| YKL045W   | PRI2      | No stall |  -0.280|
| YER055C   | HIS1      | No stall |  -0.280|
| YJL130C   | URA2      | No stall |  -0.280|
| YGL130W   | CEG1      | No stall |  -0.280|
| YGR061C   | ADE6      | No stall |  -0.281|
| YOR106W   | VAM3      | No stall |  -0.281|
| YMR199W   | CLN1      | No stall |  -0.281|
| YJL180C   | ATP12     | No stall |  -0.281|
| YLL029W   | FRA1      | Stall    |  -0.282|
| YPL244C   | HUT1      | No stall |  -0.282|
| YPR166C   | MRP2      | No stall |  -0.282|
| YDR219C   | MFB1      | No stall |  -0.282|
| YKR029C   | SET3      | Stall    |  -0.283|
| YBR126C   | TPS1      | No stall |  -0.283|
| YNL297C   | MON2      | No stall |  -0.284|
| YFR016C   | NA        | Stall    |  -0.284|
| YMR212C   | EFR3      | No stall |  -0.284|
| YGL206C   | CHC1      | No stall |  -0.285|
| YOR016C   | ERP4      | No stall |  -0.285|
| YPL135W   | ISU1      | No stall |  -0.285|
| YMR165C   | PAH1      | No stall |  -0.285|
| YBR103W   | SIF2      | No stall |  -0.285|
| YAL056W   | GPB2      | No stall |  -0.285|
| YBL022C   | PIM1      | No stall |  -0.285|
| YKL074C   | MUD2      | No stall |  -0.286|
| YML008C   | ERG6      | No stall |  -0.286|
| YLR292C   | SEC72     | No stall |  -0.287|
| YDR166C   | SEC5      | No stall |  -0.287|
| YPL145C   | KES1      | No stall |  -0.288|
| YPR191W   | QCR2      | No stall |  -0.288|
| YDR045C   | RPC11     | No stall |  -0.288|
| YLR033W   | RSC58     | Stall    |  -0.289|
| YFL008W   | SMC1      | Stall    |  -0.290|
| YPR107C   | YTH1      | No stall |  -0.290|
| YNL315C   | ATP11     | No stall |  -0.290|
| YLR367W   | RPS22B    | No stall |  -0.291|
| YFR048W   | RMD8      | No stall |  -0.291|
| YDR527W   | RBA50     | Stall    |  -0.291|
| YER014W   | HEM14     | No stall |  -0.292|
| YFR052W   | RPN12     | No stall |  -0.292|
| YNR045W   | PET494    | No stall |  -0.292|
| YDL120W   | YFH1      | No stall |  -0.293|
| YJL209W   | CBP1      | No stall |  -0.293|
| YGL048C   | RPT6      | No stall |  -0.294|
| YKL214C   | YRA2      | No stall |  -0.294|
| YHR154W   | RTT107    | No stall |  -0.294|
| YPL063W   | TIM50     | Stall    |  -0.294|
| YBR175W   | SWD3      | No stall |  -0.294|
| YDR479C   | PEX29     | No stall |  -0.295|
| YOR153W   | PDR5      | Stall    |  -0.295|
| YLR274W   | MCM5      | No stall |  -0.295|
| YER087W   | AIM10     | No stall |  -0.295|
| YDR362C   | TFC6      | Stall    |  -0.296|
| YBR168W   | PEX32     | Stall    |  -0.296|
| YLR356W   | ATG33     | No stall |  -0.296|
| YDR138W   | HPR1      | Stall    |  -0.296|
| YLR371W   | ROM2      | No stall |  -0.296|
| YCR083W   | TRX3      | No stall |  -0.296|
| YHL004W   | MRP4      | No stall |  -0.298|
| YLR254C   | NDL1      | No stall |  -0.298|
| YLR268W   | SEC22     | No stall |  -0.298|
| YLR323C   | CWC24     | No stall |  -0.299|
| YDR416W   | SYF1      | No stall |  -0.299|
| YGR287C   | IMA1      | No stall |  -0.300|
| YKL160W   | ELF1      | Stall    |  -0.301|
| YML021C   | UNG1      | No stall |  -0.301|
| YER022W   | SRB4      | No stall |  -0.301|
| YAL033W   | POP5      | No stall |  -0.301|
| YGR116W   | SPT6      | Stall    |  -0.301|
| YBL037W   | APL3      | No stall |  -0.301|
| YBR125C   | PTC4      | No stall |  -0.301|
| YBR205W   | KTR3      | No stall |  -0.302|
| YOR317W   | FAA1      | No stall |  -0.302|
| YOR064C   | YNG1      | Stall    |  -0.302|
| YNL223W   | ATG4      | No stall |  -0.302|
| YDR329C   | PEX3      | No stall |  -0.302|
| YOL032W   | OPI10     | No stall |  -0.302|
| YJR032W   | CPR7      | No stall |  -0.302|
| YER050C   | RSM18     | No stall |  -0.303|
| YOR175C   | ALE1      | No stall |  -0.303|
| YGR253C   | PUP2      | No stall |  -0.303|
| YLL022C   | HIF1      | No stall |  -0.303|
| YPR004C   | AIM45     | No stall |  -0.304|
| YBR281C   | DUG2      | No stall |  -0.304|
| YJL173C   | RFA3      | No stall |  -0.304|
| YER017C   | AFG3      | No stall |  -0.304|
| YNL281W   | HCH1      | No stall |  -0.305|
| YKL117W   | SBA1      | No stall |  -0.305|
| YLR025W   | SNF7      | No stall |  -0.306|
| YKL201C   | MNN4      | Stall    |  -0.306|
| YKL010C   | UFD4      | Stall    |  -0.306|
| YNL085W   | MKT1      | No stall |  -0.307|
| YPR103W   | PRE2      | No stall |  -0.307|
| YOR262W   | GPN2      | No stall |  -0.307|
| YOR137C   | SIA1      | No stall |  -0.307|
| YHR104W   | GRE3      | No stall |  -0.307|
| YMR260C   | TIF11     | Stall    |  -0.307|
| YBR177C   | EHT1      | No stall |  -0.307|
| YIR021W   | MRS1      | No stall |  -0.307|
| YPR120C   | CLB5      | No stall |  -0.308|
| YNL056W   | OCA2      | No stall |  -0.309|
| YKL080W   | VMA5      | No stall |  -0.309|
| YOR132W   | VPS17     | No stall |  -0.310|
| YDL191W   | RPL35A    | Stall    |  -0.310|
| YBL076C   | ILS1      | Stall    |  -0.310|
| YOR276W   | CAF20     | No stall |  -0.311|
| YOR250C   | CLP1      | No stall |  -0.311|
| YER027C   | GAL83     | No stall |  -0.312|
| YFL036W   | RPO41     | No stall |  -0.312|
| YBR052C   | RFS1      | No stall |  -0.312|
| YLR368W   | MDM30     | No stall |  -0.313|
| YJR136C   | TTI2      | No stall |  -0.315|
| YKR064W   | OAF3      | Stall    |  -0.315|
| YMR168C   | CEP3      | No stall |  -0.315|
| YML049C   | RSE1      | No stall |  -0.315|
| YDL066W   | IDP1      | No stall |  -0.316|
| YDL058W   | USO1      | No stall |  -0.316|
| YHR207C   | SET5      | No stall |  -0.316|
| YPR135W   | CTF4      | No stall |  -0.317|
| YLR148W   | PEP3      | No stall |  -0.317|
| YLR390W   | ECM19     | Stall    |  -0.318|
| YER113C   | TMN3      | No stall |  -0.318|
| YPL009C   | RQC2      | Stall    |  -0.318|
| YOR052C   | TMC1      | Stall    |  -0.319|
| YLR406C   | RPL31B    | Stall    |  -0.320|
| YEL027W   | VMA3      | No stall |  -0.320|
| YKR068C   | BET3      | No stall |  -0.321|
| YGL094C   | PAN2      | No stall |  -0.321|
| YJR094W-A | RPL43B    | No stall |  -0.321|
| YOR017W   | PET127    | No stall |  -0.321|
| YOR075W   | UFE1      | No stall |  -0.322|
| YBR208C   | DUR1,2    | Stall    |  -0.323|
| YDR244W   | PEX5      | No stall |  -0.323|
| YNL232W   | CSL4      | No stall |  -0.324|
| YKL064W   | MNR2      | Stall    |  -0.324|
| YKL171W   | NNK1      | No stall |  -0.324|
| YBR073W   | RDH54     | No stall |  -0.324|
| YDL154W   | MSH5      | No stall |  -0.324|
| YDR284C   | DPP1      | No stall |  -0.325|
| YKL016C   | ATP7      | No stall |  -0.325|
| YBR014C   | GRX7      | No stall |  -0.326|
| YDR517W   | GRH1      | Stall    |  -0.326|
| YKL104C   | GFA1      | No stall |  -0.326|
| YNL107W   | YAF9      | No stall |  -0.327|
| YLR220W   | CCC1      | No stall |  -0.327|
| YPR165W   | RHO1      | Stall    |  -0.327|
| YDR407C   | TRS120    | No stall |  -0.327|
| YBR189W   | RPS9B     | No stall |  -0.327|
| YNL286W   | CUS2      | No stall |  -0.328|
| YDR515W   | SLF1      | Stall    |  -0.328|
| YCL035C   | GRX1      | No stall |  -0.328|
| YDL164C   | CDC9      | No stall |  -0.329|
| YGL006W   | PMC1      | Stall    |  -0.330|
| YMR092C   | AIP1      | No stall |  -0.330|
| YMR231W   | PEP5      | No stall |  -0.330|
| YGR266W   | NA        | No stall |  -0.330|
| YOR298C-A | MBF1      | No stall |  -0.331|
| YJL146W   | IDS2      | No stall |  -0.332|
| YMR163C   | INP2      | No stall |  -0.332|
| YHR017W   | YSC83     | No stall |  -0.332|
| YBL047C   | EDE1      | No stall |  -0.332|
| YJL047C   | RTT101    | No stall |  -0.334|
| YGR048W   | UFD1      | Stall    |  -0.335|
| YML011C   | RAD33     | No stall |  -0.335|
| YMR195W   | ICY1      | No stall |  -0.335|
| YDL146W   | LDB17     | Stall    |  -0.336|
| YOR110W   | TFC7      | No stall |  -0.336|
| YBR115C   | LYS2      | No stall |  -0.336|
| YOL127W   | RPL25     | No stall |  -0.336|
| YJR122W   | IBA57     | No stall |  -0.336|
| YHR167W   | THP2      | No stall |  -0.336|
| YJL005W   | CYR1      | Stall    |  -0.337|
| YKL092C   | BUD2      | Stall    |  -0.337|
| YPR073C   | LTP1      | No stall |  -0.337|
| YMR112C   | MED11     | No stall |  -0.337|
| YGL093W   | SPC105    | No stall |  -0.337|
| YER013W   | PRP22     | Stall    |  -0.338|
| YDR484W   | VPS52     | No stall |  -0.338|
| YHR198C   | AIM18     | No stall |  -0.338|
| YNL067W   | RPL9B     | No stall |  -0.339|
| YLR447C   | VMA6      | No stall |  -0.339|
| YKL113C   | RAD27     | No stall |  -0.339|
| YKR054C   | DYN1      | No stall |  -0.339|
| YBR165W   | UBS1      | No stall |  -0.339|
| YKL126W   | YPK1      | No stall |  -0.339|
| YOR125C   | CAT5      | No stall |  -0.340|
| YDR350C   | ATP22     | No stall |  -0.340|
| YHR203C   | RPS4A     | No stall |  -0.340|
| YER074W-A | YOS1      | No stall |  -0.341|
| YPR134W   | MSS18     | No stall |  -0.341|
| YPR060C   | ARO7      | No stall |  -0.342|
| YBR009C   | HHF2      | No stall |  -0.342|
| YBR150C   | TBS1      | Stall    |  -0.342|
| YPL249C-A | RPL36B    | No stall |  -0.343|
| YKL085W   | MDH1      | No stall |  -0.343|
| YNL151C   | RPC31     | No stall |  -0.343|
| YJL057C   | IKS1      | Stall    |  -0.343|
| YHL010C   | ETP1      | Stall    |  -0.344|
| YBR025C   | OLA1      | No stall |  -0.344|
| YPL161C   | BEM4      | No stall |  -0.345|
| YDR529C   | QCR7      | No stall |  -0.346|
| YDL126C   | CDC48     | No stall |  -0.346|
| YDL018C   | ERP3      | Stall    |  -0.347|
| YML036W   | CGI121    | No stall |  -0.347|
| YPL178W   | CBC2      | No stall |  -0.347|
| YGL157W   | ARI1      | No stall |  -0.347|
| YOL123W   | HRP1      | No stall |  -0.347|
| YOR165W   | SEY1      | No stall |  -0.348|
| YNL224C   | SQS1      | Stall    |  -0.348|
| YBL003C   | HTA2      | No stall |  -0.348|
| YJL124C   | LSM1      | No stall |  -0.348|
| YGL241W   | KAP114    | No stall |  -0.349|
| YNL252C   | MRPL17    | No stall |  -0.349|
| YKL002W   | DID4      | No stall |  -0.350|
| YER162C   | RAD4      | Stall    |  -0.351|
| YLR330W   | CHS5      | Stall    |  -0.351|
| YNL041C   | COG6      | No stall |  -0.351|
| YCR014C   | POL4      | No stall |  -0.351|
| YDR531W   | CAB1      | No stall |  -0.351|
| YHR121W   | LSM12     | No stall |  -0.352|
| YLR272C   | YCS4      | No stall |  -0.352|
| YOR280C   | FSH3      | No stall |  -0.352|
| YOR368W   | RAD17     | No stall |  -0.353|
| YGR012W   | MCY1      | No stall |  -0.353|
| YOR254C   | SEC63     | No stall |  -0.353|
| YLR077W   | FMP25     | No stall |  -0.354|
| YDR452W   | PPN1      | Stall    |  -0.355|
| YNL257C   | SIP3      | No stall |  -0.356|
| YKL140W   | TGL1      | No stall |  -0.356|
| YPL206C   | PGC1      | No stall |  -0.356|
| YDL227C   | HO        | No stall |  -0.357|
| YPR149W   | NCE102    | No stall |  -0.357|
| YJL031C   | BET4      | No stall |  -0.357|
| YMR179W   | SPT21     | No stall |  -0.358|
| YNL084C   | END3      | Stall    |  -0.358|
| YOR164C   | GET4      | No stall |  -0.358|
| YDR063W   | AIM7      | No stall |  -0.359|
| YJR145C   | RPS4A     | No stall |  -0.359|
| YDR100W   | TVP15     | No stall |  -0.359|
| YJL157C   | FAR1      | No stall |  -0.359|
| YIL070C   | MAM33     | No stall |  -0.360|
| YMR023C   | MSS1      | No stall |  -0.360|
| YKR087C   | OMA1      | No stall |  -0.361|
| YNL166C   | BNI5      | No stall |  -0.362|
| YJL071W   | ARG2      | Stall    |  -0.363|
| YML034W   | SRC1      | Stall    |  -0.363|
| YLR383W   | SMC6      | No stall |  -0.364|
| YDR320C-A | DAD4      | No stall |  -0.364|
| YCR082W   | AHC2      | No stall |  -0.364|
| YNL077W   | APJ1      | Stall    |  -0.365|
| YOR101W   | RAS1      | No stall |  -0.365|
| YMR267W   | PPA2      | No stall |  -0.366|
| YJL139C   | YUR1      | No stall |  -0.366|
| YLR200W   | YKE2      | No stall |  -0.367|
| YJL099W   | CHS6      | No stall |  -0.367|
| YPL120W   | VPS30     | No stall |  -0.368|
| YPR056W   | TFB4      | Stall    |  -0.368|
| YJR055W   | HIT1      | No stall |  -0.369|
| YGL027C   | CWH41     | No stall |  -0.369|
| YDR229W   | IVY1      | No stall |  -0.369|
| YKL179C   | COY1      | No stall |  -0.369|
| YKL114C   | APN1      | Stall    |  -0.370|
| YOR104W   | PIN2      | Stall    |  -0.370|
| YLR066W   | SPC3      | No stall |  -0.371|
| YHR090C   | YNG2      | Stall    |  -0.371|
| YER174C   | GRX4      | No stall |  -0.371|
| YOR039W   | CKB2      | No stall |  -0.371|
| YPL050C   | MNN9      | No stall |  -0.373|
| YGL247W   | BRR6      | No stall |  -0.374|
| YFR015C   | GSY1      | No stall |  -0.374|
| YPL210C   | SRP72     | Stall    |  -0.374|
| YIR004W   | DJP1      | No stall |  -0.375|
| YDR469W   | SDC1      | No stall |  -0.375|
| YKR084C   | HBS1      | Stall    |  -0.376|
| YPL094C   | SEC62     | Stall    |  -0.376|
| YMR044W   | IOC4      | Stall    |  -0.377|
| YCR044C   | PER1      | No stall |  -0.377|
| YHR141C   | RPL42A    | No stall |  -0.377|
| YAL055W   | PEX22     | Stall    |  -0.377|
| YJL167W   | ERG20     | No stall |  -0.377|
| YPR036W-A | SPO24     | No stall |  -0.377|
| YMR072W   | ABF2      | No stall |  -0.378|
| YFR010W   | UBP6      | Stall    |  -0.378|
| YJL102W   | MEF2      | No stall |  -0.380|
| YGR094W   | VAS1      | Stall    |  -0.380|
| YBR248C   | HIS7      | No stall |  -0.380|
| YGL226W   | MTC3      | No stall |  -0.381|
| YJR092W   | BUD4      | No stall |  -0.381|
| YMR142C   | RPL13B    | Stall    |  -0.382|
| YLR197W   | NOP56     | Stall    |  -0.382|
| YNL048W   | ALG11     | No stall |  -0.382|
| YML067C   | ERV41     | No stall |  -0.383|
| YHR068W   | DYS1      | No stall |  -0.383|
| YOR044W   | IRC23     | No stall |  -0.383|
| YOR193W   | PEX27     | No stall |  -0.383|
| YJL140W   | RPB4      | Stall    |  -0.384|
| YLR259C   | HSP60     | No stall |  -0.384|
| YLR045C   | STU2      | No stall |  -0.384|
| YDR167W   | TAF10     | No stall |  -0.385|
| YIL116W   | HIS5      | No stall |  -0.385|
| YPR102C   | RPL11A    | Stall    |  -0.387|
| YIL106W   | MOB1      | No stall |  -0.388|
| YLR351C   | NIT3      | No stall |  -0.388|
| YLR424W   | SPP382    | Stall    |  -0.388|
| YLR248W   | RCK2      | No stall |  -0.388|
| YMR251W-A | HOR7      | No stall |  -0.388|
| YOL149W   | DCP1      | No stall |  -0.389|
| YFL028C   | CAF16     | No stall |  -0.389|
| YIL078W   | THS1      | No stall |  -0.390|
| YPR029C   | APL4      | No stall |  -0.390|
| YKL193C   | SDS22     | No stall |  -0.390|
| YDR495C   | VPS3      | No stall |  -0.391|
| YDR051C   | DET1      | No stall |  -0.391|
| YHR016C   | YSC84     | Stall    |  -0.391|
| YIL146C   | ATG32     | Stall    |  -0.392|
| YDR502C   | SAM2      | No stall |  -0.392|
| YGL058W   | RAD6      | Stall    |  -0.392|
| YMR136W   | GAT2      | Stall    |  -0.393|
| YAL009W   | SPO7      | No stall |  -0.394|
| YLR107W   | REX3      | No stall |  -0.394|
| YLR396C   | VPS33     | No stall |  -0.394|
| YDL045W-A | MRP10     | No stall |  -0.395|
| YKR016W   | MIC60     | No stall |  -0.395|
| YKL027W   | TCD2      | Stall    |  -0.396|
| YGL229C   | SAP4      | Stall    |  -0.396|
| YDR072C   | IPT1      | Stall    |  -0.397|
| YHR004C   | NEM1      | No stall |  -0.397|
| YDL178W   | DLD2      | No stall |  -0.397|
| YNL081C   | SWS2      | No stall |  -0.398|
| YIL090W   | ICE2      | No stall |  -0.398|
| YOR258W   | HNT3      | No stall |  -0.398|
| YDL171C   | GLT1      | No stall |  -0.400|
| YPL139C   | UME1      | No stall |  -0.400|
| YGR275W   | RTT102    | No stall |  -0.400|
| YDL098C   | SNU23     | No stall |  -0.400|
| YLR309C   | IMH1      | Stall    |  -0.400|
| YNL330C   | RPD3      | No stall |  -0.400|
| YNL206C   | RTT106    | No stall |  -0.401|
| YKL091C   | NA        | No stall |  -0.401|
| YOR321W   | PMT3      | No stall |  -0.401|
| YMR304W   | UBP15     | No stall |  -0.401|
| YLR191W   | PEX13     | No stall |  -0.402|
| YPL219W   | PCL8      | No stall |  -0.402|
| YNL088W   | TOP2      | Stall    |  -0.402|
| YDR294C   | DPL1      | No stall |  -0.402|
| YJR144W   | MGM101    | Stall    |  -0.402|
| YHR023W   | MYO1      | No stall |  -0.402|
| YNR022C   | MRPL50    | No stall |  -0.403|
| YMR005W   | TAF4      | Stall    |  -0.404|
| YBR230C   | OM14      | No stall |  -0.404|
| YPL153C   | RAD53     | No stall |  -0.405|
| YGL233W   | SEC15     | No stall |  -0.405|
| YFR033C   | QCR6      | No stall |  -0.405|
| YGL019W   | CKB1      | No stall |  -0.405|
| YOL069W   | NUF2      | No stall |  -0.405|
| YBL087C   | RPL23A    | Stall    |  -0.405|
| YGR150C   | CCM1      | No stall |  -0.405|
| YHR117W   | TOM71     | Stall    |  -0.406|
| YJL126W   | NIT2      | No stall |  -0.406|
| YGL161C   | YIP5      | No stall |  -0.406|
| YBR162W-A | YSY6      | Stall    |  -0.406|
| YMR152W   | YIM1      | No stall |  -0.407|
| YOR076C   | SKI7      | No stall |  -0.407|
| YJL074C   | SMC3      | No stall |  -0.408|
| YMR001C   | CDC5      | No stall |  -0.408|
| YPL265W   | DIP5      | No stall |  -0.409|
| YMR089C   | YTA12     | Stall    |  -0.410|
| YKR057W   | RPS21A    | No stall |  -0.410|
| YPL218W   | SAR1      | No stall |  -0.410|
| YDL083C   | RPS16A    | No stall |  -0.411|
| YLR210W   | CLB4      | Stall    |  -0.411|
| YHL038C   | CBP2      | No stall |  -0.411|
| YEL029C   | BUD16     | No stall |  -0.412|
| YDR146C   | SWI5      | Stall    |  -0.412|
| YPL104W   | MSD1      | No stall |  -0.412|
| YPL158C   | AIM44     | Stall    |  -0.413|
| YBR185C   | MBA1      | No stall |  -0.413|
| YJL085W   | EXO70     | No stall |  -0.413|
| YBR039W   | ATP3      | No stall |  -0.414|
| YPL003W   | ULA1      | No stall |  -0.414|
| YML109W   | ZDS2      | Stall    |  -0.414|
| YML077W   | BET5      | No stall |  -0.415|
| YDR022C   | ATG31     | No stall |  -0.416|
| YMR198W   | CIK1      | Stall    |  -0.416|
| YHR144C   | DCD1      | No stall |  -0.416|
| YJL112W   | MDV1      | Stall    |  -0.417|
| YMR083W   | ADH3      | No stall |  -0.417|
| YKL007W   | CAP1      | No stall |  -0.418|
| YDR258C   | HSP78     | No stall |  -0.418|
| YDR194C   | MSS116    | No stall |  -0.418|
| YGR109C   | CLB6      | No stall |  -0.419|
| YMR004W   | MVP1      | No stall |  -0.419|
| YPL231W   | FAS2      | No stall |  -0.419|
| YIL158W   | AIM20     | No stall |  -0.420|
| YDR068W   | DOS2      | No stall |  -0.421|
| YLR382C   | NAM2      | No stall |  -0.421|
| YOR087W   | YVC1      | Stall    |  -0.422|
| YPR139C   | LOA1      | Stall    |  -0.424|
| YGR255C   | COQ6      | No stall |  -0.424|
| YMR191W   | SPG5      | No stall |  -0.424|
| YFR003C   | YPI1      | Stall    |  -0.424|
| YKR098C   | UBP11     | Stall    |  -0.425|
| YBR218C   | PYC2      | Stall    |  -0.425|
| YKL141W   | SDH3      | No stall |  -0.426|
| YML078W   | CPR3      | No stall |  -0.426|
| YDR518W   | EUG1      | No stall |  -0.426|
| YDL029W   | ARP2      | No stall |  -0.426|
| YLR441C   | RPS1A     | Stall    |  -0.426|
| YMR060C   | SAM37     | No stall |  -0.426|
| YDL232W   | OST4      | No stall |  -0.427|
| YJR074W   | MOG1      | No stall |  -0.427|
| YKL148C   | SDH1      | No stall |  -0.427|
| YLR442C   | SIR3      | Stall    |  -0.427|
| YOR279C   | RFM1      | Stall    |  -0.427|
| YOL042W   | NGL1      | No stall |  -0.429|
| YOR211C   | MGM1      | No stall |  -0.430|
| YBR088C   | POL30     | No stall |  -0.430|
| YOR367W   | SCP1      | No stall |  -0.431|
| YGL037C   | PNC1      | No stall |  -0.431|
| YDR071C   | PAA1      | No stall |  -0.431|
| YOR018W   | ROD1      | No stall |  -0.432|
| YOL142W   | RRP40     | No stall |  -0.432|
| YKL167C   | MRP49     | No stall |  -0.432|
| YDR345C   | HXT3      | No stall |  -0.433|
| YJR050W   | ISY1      | Stall    |  -0.433|
| YFR031C   | SMC2      | No stall |  -0.433|
| YNR030W   | ALG12     | No stall |  -0.433|
| YLR448W   | RPL6B     | Stall    |  -0.434|
| YKL065C   | YET1      | No stall |  -0.434|
| YDL075W   | RPL31A    | Stall    |  -0.434|
| YHR057C   | CPR2      | No stall |  -0.435|
| YMR024W   | MRPL3     | No stall |  -0.435|
| YDR388W   | RVS167    | No stall |  -0.435|
| YDR375C   | BCS1      | No stall |  -0.435|
| YDR005C   | MAF1      | No stall |  -0.436|
| YDR393W   | SHE9      | No stall |  -0.436|
| YPL259C   | APM1      | Stall    |  -0.437|
| YER117W   | RPL23A    | Stall    |  -0.437|
| YHR028C   | DAP2      | No stall |  -0.437|
| YNL142W   | MEP2      | No stall |  -0.438|
| YDL107W   | MSS2      | No stall |  -0.438|
| YIL009W   | FAA3      | Stall    |  -0.438|
| YGR180C   | RNR4      | No stall |  -0.439|
| YDR382W   | RPP2B     | No stall |  -0.439|
| YDR289C   | RTT103    | No stall |  -0.441|
| YJL051W   | IRC8      | No stall |  -0.442|
| YIL045W   | PIG2      | Stall    |  -0.442|
| YLR420W   | URA4      | No stall |  -0.442|
| YPR124W   | CTR1      | No stall |  -0.442|
| YPR128C   | ANT1      | No stall |  -0.443|
| YER141W   | COX15     | No stall |  -0.443|
| YEL050C   | RML2      | Stall    |  -0.444|
| YGL055W   | OLE1      | Stall    |  -0.444|
| YBR161W   | CSH1      | No stall |  -0.445|
| YDR204W   | COQ4      | Stall    |  -0.445|
| YLR262C-A | TMA7      | Stall    |  -0.446|
| YHR001W-A | QCR10     | No stall |  -0.447|
| YPL100W   | ATG21     | No stall |  -0.447|
| YNL289W   | PCL1      | No stall |  -0.447|
| YDR525W-A | SNA2      | No stall |  -0.447|
| YDL004W   | ATP16     | No stall |  -0.448|
| YNR007C   | ATG3      | Stall    |  -0.448|
| YEL051W   | VMA8      | No stall |  -0.450|
| YIL105C   | SLM1      | No stall |  -0.450|
| YLR189C   | ATG26     | Stall    |  -0.451|
| YOR096W   | RPS7A     | No stall |  -0.452|
| YLR051C   | FCF2      | Stall    |  -0.453|
| YFL025C   | BST1      | No stall |  -0.454|
| YBR256C   | RIB5      | No stall |  -0.455|
| YMR115W   | MGR3      | No stall |  -0.455|
| YOL090W   | MSH2      | Stall    |  -0.456|
| YKL206C   | ADD66     | No stall |  -0.457|
| YBR058C-A | TSC3      | No stall |  -0.457|
| Q0085     | ATP6      | No stall |  -0.457|
| YOL034W   | SMC5      | No stall |  -0.457|
| YPL111W   | CAR1      | No stall |  -0.457|
| YLR240W   | VPS34     | No stall |  -0.459|
| YCR019W   | MAK32     | No stall |  -0.460|
| YOR265W   | RBL2      | No stall |  -0.460|
| YNL231C   | PDR16     | No stall |  -0.461|
| YER051W   | JHD1      | Stall    |  -0.464|
| YPL064C   | CWC27     | No stall |  -0.464|
| YBR221C   | PDB1      | No stall |  -0.464|
| YOR117W   | RPT5      | No stall |  -0.465|
| YMR287C   | DSS1      | Stall    |  -0.466|
| YOL109W   | ZEO1      | No stall |  -0.467|
| YDR397C   | NCB2      | No stall |  -0.468|
| YGL054C   | ERV14     | No stall |  -0.468|
| YHL030W   | ECM29     | No stall |  -0.470|
| YBR139W   | NA        | Stall    |  -0.471|
| YPL010W   | RET3      | No stall |  -0.471|
| YLL040C   | VPS13     | Stall    |  -0.472|
| YOR196C   | LIP5      | Stall    |  -0.473|
| YNL026W   | SAM50     | No stall |  -0.474|
| YDL091C   | UBX3      | No stall |  -0.474|
| YOR142W   | LSC1      | No stall |  -0.474|
| YLR372W   | ELO3      | Stall    |  -0.475|
| YBL050W   | SEC17     | No stall |  -0.475|
| YHR001W   | OSH7      | No stall |  -0.476|
| YKL218C   | SRY1      | No stall |  -0.476|
| YGR030C   | POP6      | No stall |  -0.476|
| YBR077C   | SLM4      | No stall |  -0.477|
| YNL008C   | ASI3      | Stall    |  -0.478|
| YOR115C   | TRS33     | No stall |  -0.478|
| YGR055W   | MUP1      | No stall |  -0.478|
| YPL081W   | RPS9A     | No stall |  -0.479|
| YPR132W   | RPS23B    | No stall |  -0.480|
| YDR356W   | SPC110    | No stall |  -0.481|
| YKL056C   | TMA19     | No stall |  -0.481|
| YMR121C   | RPL15B    | Stall    |  -0.481|
| YOR357C   | SNX3      | Stall    |  -0.481|
| YNL162W   | RPL42A    | No stall |  -0.481|
| YCL024W   | KCC4      | Stall    |  -0.481|
| YOR257W   | CDC31     | No stall |  -0.483|
| YLR285W   | NNT1      | No stall |  -0.483|
| YNL133C   | FYV6      | Stall    |  -0.483|
| YMR311C   | GLC8      | No stall |  -0.484|
| YDR447C   | RPS17B    | Stall    |  -0.484|
| YPR083W   | MDM36     | Stall    |  -0.484|
| YBL036C   | NA        | No stall |  -0.484|
| YJL155C   | FBP26     | Stall    |  -0.485|
| YEL022W   | GEA2      | No stall |  -0.485|
| YFL018C   | LPD1      | No stall |  -0.485|
| YBR132C   | AGP2      | No stall |  -0.485|
| YOR136W   | IDH2      | No stall |  -0.485|
| YAL030W   | SNC1      | No stall |  -0.486|
| YOR143C   | THI80     | No stall |  -0.488|
| YKL190W   | CNB1      | No stall |  -0.488|
| YDR507C   | GIN4      | Stall    |  -0.489|
| YIL017C   | VID28     | No stall |  -0.489|
| YBL030C   | PET9      | No stall |  -0.489|
| YOR122C   | PFY1      | No stall |  -0.489|
| YNL258C   | DSL1      | No stall |  -0.489|
| YPL271W   | ATP15     | No stall |  -0.490|
| YOL065C   | INP54     | No stall |  -0.490|
| YLL041C   | SDH2      | No stall |  -0.490|
| YCR036W   | RBK1      | No stall |  -0.490|
| YPR006C   | ICL2      | No stall |  -0.491|
| YNL121C   | TOM70     | No stall |  -0.491|
| YKR085C   | MRPL20    | Stall    |  -0.491|
| YOL064C   | MET22     | No stall |  -0.493|
| YBR057C   | MUM2      | No stall |  -0.493|
| YMR257C   | PET111    | No stall |  -0.493|
| YOL086C   | ADH1      | No stall |  -0.493|
| YFL004W   | VTC2      | Stall    |  -0.494|
| YIL004C   | BET1      | No stall |  -0.494|
| YOL023W   | IFM1      | No stall |  -0.494|
| YER023W   | PRO3      | No stall |  -0.495|
| YOR147W   | MDM32     | No stall |  -0.496|
| YOR025W   | HST3      | No stall |  -0.498|
| YMR150C   | IMP1      | No stall |  -0.499|
| YLR103C   | CDC45     | Stall    |  -0.500|
| YBR084C-A | RPL19B    | No stall |  -0.500|
| YDR308C   | SRB7      | No stall |  -0.501|
| YLR105C   | SEN2      | No stall |  -0.501|
| YOR354C   | MSC6      | No stall |  -0.502|
| YBR129C   | OPY1      | Stall    |  -0.502|
| YKR101W   | SIR1      | Stall    |  -0.503|
| YDR378C   | LSM6      | No stall |  -0.503|
| YIL125W   | KGD1      | No stall |  -0.504|
| YMR241W   | YHM2      | No stall |  -0.504|
| YDR238C   | SEC26     | Stall    |  -0.504|
| YPR148C   | NA        | No stall |  -0.505|
| YDR118W   | APC4      | No stall |  -0.505|
| YDR322C-A | TIM11     | No stall |  -0.505|
| YIL026C   | IRR1      | Stall    |  -0.507|
| YJR103W   | URA8      | No stall |  -0.507|
| YMR302C   | YME2      | No stall |  -0.507|
| YOR212W   | STE4      | No stall |  -0.508|
| YNL287W   | SEC21     | No stall |  -0.508|
| YDR516C   | EMI2      | No stall |  -0.509|
| YNL005C   | MRP7      | No stall |  -0.510|
| YDL081C   | RPP1A     | No stall |  -0.511|
| YJL044C   | GYP6      | No stall |  -0.512|
| YBR037C   | SCO1      | No stall |  -0.512|
| YBR043C   | QDR3      | No stall |  -0.513|
| YIL068C   | SEC6      | No stall |  -0.514|
| YCR024C-A | PMP1      | No stall |  -0.514|
| YBR127C   | VMA2      | No stall |  -0.515|
| YMR104C   | YPK2      | No stall |  -0.517|
| YJL003W   | COX16     | No stall |  -0.520|
| YGR282C   | BGL2      | No stall |  -0.520|
| YBR234C   | ARC40     | No stall |  -0.521|
| YGR183C   | QCR9      | No stall |  -0.522|
| YBL027W   | RPL19B    | No stall |  -0.522|
| YMR140W   | SIP5      | No stall |  -0.524|
| YEL061C   | CIN8      | No stall |  -0.525|
| YMR202W   | ERG2      | No stall |  -0.525|
| YLR239C   | LIP2      | No stall |  -0.526|
| YPR062W   | FCY1      | No stall |  -0.526|
| YER080W   | AIM9      | No stall |  -0.526|
| YDL092W   | SRP14     | Stall    |  -0.528|
| YGR086C   | PIL1      | No stall |  -0.528|
| YOR222W   | ODC2      | No stall |  -0.529|
| YNL129W   | NRK1      | Stall    |  -0.529|
| YJR048W   | CYC1      | No stall |  -0.529|
| YGL116W   | CDC20     | No stall |  -0.529|
| YPL042C   | SSN3      | No stall |  -0.530|
| YCR024C   | SLM5      | No stall |  -0.530|
| YML055W   | SPC2      | No stall |  -0.532|
| YOR037W   | CYC2      | No stall |  -0.532|
| YNL082W   | PMS1      | No stall |  -0.532|
| YGR288W   | MAL13     | No stall |  -0.533|
| YDL010W   | GRX6      | No stall |  -0.533|
| YOR259C   | RPT4      | No stall |  -0.535|
| YJR139C   | HOM6      | No stall |  -0.535|
| YOL146W   | PSF3      | No stall |  -0.538|
| YMR009W   | ADI1      | No stall |  -0.538|
| YMR173W   | DDR48     | No stall |  -0.538|
| YDR181C   | SAS4      | Stall    |  -0.539|
| YGR244C   | LSC2      | No stall |  -0.539|
| YIL124W   | AYR1      | No stall |  -0.539|
| YDL149W   | ATG9      | No stall |  -0.539|
| YDR046C   | BAP3      | No stall |  -0.541|
| YLR093C   | NYV1      | No stall |  -0.542|
| YBR199W   | KTR4      | No stall |  -0.544|
| YOR086C   | TCB1      | Stall    |  -0.544|
| YDL003W   | MCD1      | No stall |  -0.547|
| YPR111W   | DBF20     | Stall    |  -0.548|
| YLR181C   | VTA1      | No stall |  -0.548|
| YBR211C   | AME1      | No stall |  -0.549|
| YHL028W   | WSC4      | No stall |  -0.550|
| YBR001C   | NTH2      | No stall |  -0.550|
| YCR012W   | PGK1      | No stall |  -0.551|
| YBR171W   | SEC66     | No stall |  -0.551|
| YKR035W-A | DID2      | No stall |  -0.552|
| YGL209W   | MIG2      | No stall |  -0.554|
| YPR051W   | MAK3      | No stall |  -0.555|
| YJL060W   | BNA3      | No stall |  -0.556|
| YPL078C   | ATP4      | No stall |  -0.556|
| YOL012C   | HTZ1      | No stall |  -0.556|
| YDL061C   | RPS29B    | No stall |  -0.557|
| YDR305C   | HNT2      | No stall |  -0.557|
| YIL107C   | PFK26     | Stall    |  -0.557|
| YJR110W   | YMR1      | No stall |  -0.557|
| YLR167W   | RPS31     | Stall    |  -0.558|
| YPL195W   | APL5      | Stall    |  -0.560|
| YNL213C   | RRG9      | Stall    |  -0.560|
| YOL096C   | COQ3      | No stall |  -0.560|
| YFL041W   | FET5      | No stall |  -0.561|
| YNL137C   | NAM9      | No stall |  -0.561|
| YBR131W   | CCZ1      | No stall |  -0.561|
| YER142C   | MAG1      | Stall    |  -0.561|
| YOR244W   | ESA1      | Stall    |  -0.562|
| YCL001W   | RER1      | No stall |  -0.562|
| YGR092W   | DBF2      | Stall    |  -0.562|
| YDR168W   | CDC37     | No stall |  -0.562|
| YDL202W   | MRPL11    | No stall |  -0.563|
| YPL007C   | TFC8      | No stall |  -0.564|
| YDR405W   | MRP20     | No stall |  -0.565|
| YMR097C   | MTG1      | No stall |  -0.566|
| YKL195W   | MIA40     | No stall |  -0.566|
| YGR121C   | MEP1      | No stall |  -0.567|
| YNL233W   | BNI4      | Stall    |  -0.568|
| YNR041C   | COQ2      | No stall |  -0.568|
| YDL089W   | NUR1      | Stall    |  -0.569|
| YNL126W   | SPC98     | No stall |  -0.570|
| YHR172W   | SPC97     | No stall |  -0.570|
| YKL161C   | KDX1      | No stall |  -0.570|
| YMR077C   | VPS20     | No stall |  -0.572|
| YKR067W   | GPT2      | No stall |  -0.574|
| YJR088C   | EMC2      | No stall |  -0.575|
| YGL021W   | ALK1      | No stall |  -0.576|
| YMR284W   | YKU70     | Stall    |  -0.577|
| YOR195W   | SLK19     | Stall    |  -0.578|
| YBR268W   | MRPL37    | Stall    |  -0.578|
| YDR424C   | DYN2      | No stall |  -0.579|
| YMR062C   | ARG7      | No stall |  -0.581|
| YGR047C   | TFC4      | Stall    |  -0.582|
| YKR089C   | TGL4      | No stall |  -0.582|
| YMR061W   | RNA14     | No stall |  -0.583|
| YHL001W   | RPL14B    | Stall    |  -0.583|
| YLR260W   | LCB5      | Stall    |  -0.584|
| YLR145W   | RMP1      | Stall    |  -0.584|
| YOR130C   | ORT1      | No stall |  -0.586|
| YKL155C   | RSM22     | Stall    |  -0.586|
| YKR002W   | PAP1      | Stall    |  -0.587|
| YDR196C   | CAB5      | No stall |  -0.587|
| YNL293W   | MSB3      | Stall    |  -0.587|
| YGR220C   | MRPL9     | No stall |  -0.588|
| YGR118W   | RPS23B    | No stall |  -0.588|
| YDL147W   | RPN5      | No stall |  -0.590|
| YGR135W   | PRE9      | No stall |  -0.592|
| YLR300W   | EXG1      | No stall |  -0.592|
| YBR097W   | VPS15     | No stall |  -0.593|
| YGR019W   | UGA1      | No stall |  -0.594|
| YPL092W   | SSU1      | No stall |  -0.594|
| YBR003W   | COQ1      | No stall |  -0.594|
| YDL130W-A | STF1      | No stall |  -0.595|
| YPL116W   | HOS3      | No stall |  -0.597|
| YPL183W-A | RTC6      | No stall |  -0.597|
| YGL135W   | RPL1A     | Stall    |  -0.598|
| YNL177C   | MRPL22    | No stall |  -0.598|
| YPL069C   | BTS1      | No stall |  -0.598|
| YGL248W   | PDE1      | No stall |  -0.598|
| YGR277C   | CAB4      | No stall |  -0.601|
| YOR054C   | VHS3      | No stall |  -0.601|
| YKL163W   | PIR3      | No stall |  -0.603|
| YNL192W   | CHS1      | No stall |  -0.604|
| YDR434W   | GPI17     | No stall |  -0.605|
| YMR098C   | ATP25     | Stall    |  -0.607|
| YCR028C-A | RIM1      | No stall |  -0.607|
| YNL225C   | CNM67     | No stall |  -0.609|
| YOR124C   | UBP2      | No stall |  -0.609|
| YKR088C   | TVP38     | No stall |  -0.609|
| YJR130C   | STR2      | No stall |  -0.611|
| YBR041W   | FAT1      | No stall |  -0.611|
| YDL181W   | INH1      | No stall |  -0.615|
| YLR183C   | TOS4      | No stall |  -0.616|
| YBR181C   | RPS6A     | Stall    |  -0.616|
| YJL052W   | TDH1      | No stall |  -0.616|
| YDL084W   | SUB2      | No stall |  -0.618|
| YML120C   | NDI1      | No stall |  -0.619|
| YIL027C   | EMC5      | No stall |  -0.619|
| YPL053C   | KTR6      | No stall |  -0.620|
| YBL001C   | ECM15     | No stall |  -0.620|
| YAL017W   | PSK1      | Stall    |  -0.621|
| YDR175C   | RSM24     | Stall    |  -0.621|
| YNL173C   | MDG1      | Stall    |  -0.622|
| YMR143W   | RPS16A    | No stall |  -0.622|
| YMR076C   | PDS5      | Stall    |  -0.622|
| YKL138C   | MRPL31    | Stall    |  -0.625|
| YDR497C   | ITR1      | No stall |  -0.625|
| YHR132W-A | IGO2      | No stall |  -0.625|
| YJR135W-A | TIM8      | No stall |  -0.625|
| YMR280C   | CAT8      | Stall    |  -0.627|
| YMR315W   | NA        | No stall |  -0.628|
| YDL044C   | MTF2      | No stall |  -0.628|
| YLR312W-A | MRPL15    | No stall |  -0.628|
| YJR113C   | RSM7      | No stall |  -0.629|
| YER133W   | GLC7      | No stall |  -0.629|
| YMR013C   | SEC59     | No stall |  -0.629|
| YLR369W   | SSQ1      | No stall |  -0.633|
| YOR185C   | GSP2      | No stall |  -0.634|
| YGR099W   | TEL2      | No stall |  -0.638|
| YGR179C   | OKP1      | Stall    |  -0.638|
| YMR087W   | NA        | No stall |  -0.639|
| YML030W   | RCF1      | No stall |  -0.640|
| YGL187C   | COX4      | No stall |  -0.641|
| YLR170C   | APS1      | No stall |  -0.641|
| Q0045     | COX1      | No stall |  -0.642|
| YGR029W   | ERV1      | Stall    |  -0.643|
| YDR510W   | SMT3      | No stall |  -0.644|
| YLR185W   | RPL37A    | No stall |  -0.644|
| YGL153W   | PEX14     | No stall |  -0.644|
| YLR262C   | YPT6      | No stall |  -0.644|
| YKL040C   | NFU1      | No stall |  -0.645|
| YBR137W   | NA        | No stall |  -0.646|
| YGR207C   | CIR1      | No stall |  -0.646|
| YGL226C-A | OST5      | No stall |  -0.646|
| YHL027W   | RIM101    | No stall |  -0.648|
| YLL002W   | RTT109    | No stall |  -0.648|
| YMR194W   | RPL36A    | Stall    |  -0.649|
| YGR065C   | VHT1      | No stall |  -0.650|
| YNL306W   | MRPS18    | No stall |  -0.650|
| YCR005C   | CIT2      | No stall |  -0.651|
| YBL098W   | BNA4      | No stall |  -0.651|
| YPL026C   | SKS1      | No stall |  -0.652|
| YNL169C   | PSD1      | Stall    |  -0.652|
| YKL142W   | MRP8      | No stall |  -0.654|
| YNL284C   | MRPL10    | No stall |  -0.655|
| YOR057W   | SGT1      | Stall    |  -0.655|
| YNR032C-A | HUB1      | No stall |  -0.655|
| YDR197W   | CBS2      | No stall |  -0.656|
| YOR256C   | TRE2      | No stall |  -0.657|
| YKL194C   | MST1      | No stall |  -0.658|
| YGR101W   | PCP1      | No stall |  -0.661|
| YML073C   | RPL6A     | Stall    |  -0.662|
| YDR358W   | GGA1      | No stall |  -0.663|
| YPR141C   | KAR3      | No stall |  -0.665|
| YDL065C   | PEX19     | No stall |  -0.666|
| YFR002W   | NIC96     | No stall |  -0.666|
| YMR074C   | SDD2      | No stall |  -0.667|
| YDL200C   | MGT1      | No stall |  -0.669|
| YMR282C   | AEP2      | No stall |  -0.670|
| YML060W   | OGG1      | No stall |  -0.671|
| YIL093C   | RSM25     | No stall |  -0.671|
| YDR201W   | SPC19     | Stall    |  -0.671|
| YAL007C   | ERP2      | No stall |  -0.672|
| YML012W   | ERV25     | No stall |  -0.674|
| YDR179C   | CSN9      | No stall |  -0.676|
| YER001W   | MNN1      | No stall |  -0.677|
| YOR045W   | TOM6      | No stall |  -0.678|
| YER095W   | RAD51     | No stall |  -0.678|
| YHR147C   | MRPL6     | No stall |  -0.679|
| YJL136C   | RPS21B    | No stall |  -0.680|
| YJL013C   | MAD3      | No stall |  -0.680|
| YDR041W   | RSM10     | No stall |  -0.681|
| YAL049C   | AIM2      | No stall |  -0.682|
| YPR025C   | CCL1      | No stall |  -0.683|
| YDR027C   | VPS54     | No stall |  -0.685|
| YPR155C   | NCA2      | No stall |  -0.686|
| YJR137C   | MET5      | Stall    |  -0.686|
| YHL016C   | DUR3      | No stall |  -0.687|
| YBL102W   | SFT2      | No stall |  -0.687|
| YHR037W   | PUT2      | No stall |  -0.687|
| YBR111C   | YSA1      | No stall |  -0.688|
| YKL119C   | VPH2      | Stall    |  -0.688|
| YMR117C   | SPC24     | No stall |  -0.689|
| YOL016C   | CMK2      | No stall |  -0.692|
| YDR255C   | RMD5      | Stall    |  -0.693|
| YGR147C   | NAT2      | Stall    |  -0.695|
| YMR186W   | HSC82     | Stall    |  -0.695|
| YKL013C   | ARC19     | No stall |  -0.699|
| YLR264W   | RPS28B    | No stall |  -0.699|
| YGL129C   | RSM23     | No stall |  -0.700|
| YPR158W   | CUR1      | No stall |  -0.700|
| YLR304C   | ACO1      | No stall |  -0.700|
| YGR085C   | RPL11B    | Stall    |  -0.701|
| YHR011W   | DIA4      | No stall |  -0.704|
| YLL028W   | TPO1      | No stall |  -0.705|
| YPL097W   | MSY1      | No stall |  -0.708|
| YPL090C   | RPS6A     | Stall    |  -0.711|
| YPL274W   | SAM3      | No stall |  -0.712|
| YER012W   | PRE1      | No stall |  -0.713|
| YPR160W   | GPH1      | Stall    |  -0.714|
| YOR150W   | MRPL23    | No stall |  -0.715|
| YAL003W   | EFB1      | No stall |  -0.715|
| YIR037W   | HYR1      | No stall |  -0.717|
| YLR166C   | SEC10     | No stall |  -0.717|
| YGR194C   | XKS1      | No stall |  -0.719|
| YDR381C-A | NA        | No stall |  -0.721|
| YHR060W   | VMA22     | No stall |  -0.722|
| YKL006W   | RPL14A    | Stall    |  -0.722|
| YFL005W   | SEC4      | No stall |  -0.723|
| YGL180W   | ATG1      | No stall |  -0.728|
| YBR111W-A | SUS1      | No stall |  -0.730|
| YLR079W   | SIC1      | No stall |  -0.731|
| YLR258W   | GSY2      | No stall |  -0.734|
| YIL098C   | FMC1      | No stall |  -0.734|
| YOR158W   | PET123    | No stall |  -0.734|
| YBR122C   | MRPL36    | Stall    |  -0.734|
| YMR316W   | DIA1      | Stall    |  -0.737|
| YDR342C   | HXT7      | No stall |  -0.738|
| YLR069C   | MEF1      | No stall |  -0.740|
| YGR231C   | PHB2      | No stall |  -0.740|
| YMR139W   | RIM11     | No stall |  -0.741|
| YJL166W   | QCR8      | No stall |  -0.742|
| YJL164C   | TPK1      | No stall |  -0.744|
| YBR006W   | UGA2      | No stall |  -0.744|
| YIR022W   | SEC11     | No stall |  -0.749|
| YJL082W   | IML2      | No stall |  -0.751|
| YOL143C   | RIB4      | No stall |  -0.752|
| YOR232W   | MGE1      | No stall |  -0.754|
| YMR145C   | NDE1      | No stall |  -0.754|
| YDL072C   | YET3      | No stall |  -0.755|
| Q0105     | COB       | No stall |  -0.757|
| YJL006C   | CTK2      | No stall |  -0.761|
| YOR236W   | DFR1      | No stall |  -0.761|
| YML026C   | RPS18B    | No stall |  -0.762|
| YBL090W   | MRP21     | No stall |  -0.763|
| YDR279W   | RNH202    | No stall |  -0.764|
| YEL017C-A | PMP2      | No stall |  -0.766|
| YPR002W   | PDH1      | No stall |  -0.767|
| YBL057C   | PTH2      | No stall |  -0.768|
| YDR519W   | FPR2      | No stall |  -0.769|
| YIL097W   | FYV10     | No stall |  -0.772|
| YOL136C   | PFK27     | No stall |  -0.774|
| YMR071C   | TVP18     | No stall |  -0.779|
| YMR053C   | STB2      | No stall |  -0.782|
| YDR337W   | MRPS28    | Stall    |  -0.783|
| YPL096W   | PNG1      | No stall |  -0.784|
| YER074W   | RPS24A    | Stall    |  -0.786|
| YCL030C   | HIS4      | No stall |  -0.787|
| YPR082C   | DIB1      | No stall |  -0.787|
| YDR079W   | PET100    | Stall    |  -0.791|
| YOR273C   | TPO4      | Stall    |  -0.792|
| YML110C   | COQ5      | No stall |  -0.792|
| YMR029C   | FAR8      | No stall |  -0.793|
| YDR450W   | RPS18B    | No stall |  -0.798|
| YKL124W   | SSH4      | Stall    |  -0.807|
| YBR282W   | MRPL27    | No stall |  -0.807|
| YNL265C   | IST1      | Stall    |  -0.808|
| YJR133W   | XPT1      | No stall |  -0.808|
| YPL188W   | POS5      | No stall |  -0.809|
| YGL220W   | BOL2      | No stall |  -0.809|
| YCR071C   | IMG2      | No stall |  -0.810|
| YNL030W   | HHF2      | No stall |  -0.813|
| YOR266W   | PNT1      | No stall |  -0.817|
| YDR460W   | TFB3      | Stall    |  -0.821|
| YGR027C   | RPS25A    | Stall    |  -0.823|
| YGR132C   | PHB1      | No stall |  -0.825|
| YOR187W   | TUF1      | No stall |  -0.825|
| YDR173C   | ARG82     | No stall |  -0.826|
| YDL192W   | ARF1      | No stall |  -0.829|
| YGR020C   | VMA7      | No stall |  -0.832|
| YOL053W   | AIM39     | No stall |  -0.832|
| YDR383C   | NKP1      | No stall |  -0.833|
| YHR002W   | LEU5      | No stall |  -0.835|
| YML121W   | GTR1      | No stall |  -0.837|
| YNL302C   | RPS19B    | No stall |  -0.837|
| YBL038W   | MRPL16    | No stall |  -0.843|
| YGL115W   | SNF4      | No stall |  -0.843|
| YDR084C   | TVP23     | No stall |  -0.849|
| YDR086C   | SSS1      | No stall |  -0.850|
| YDL130W   | RPP1B     | No stall |  -0.854|
| YDR494W   | RSM28     | No stall |  -0.858|
| YGR148C   | RPL24B    | Stall    |  -0.859|
| YGL059W   | PKP2      | No stall |  -0.861|
| YGR165W   | MRPS35    | No stall |  -0.862|
| YBR066C   | NRG2      | Stall    |  -0.862|
| YPR098C   | NA        | No stall |  -0.864|
| YBR146W   | MRPS9     | Stall    |  -0.864|
| YCR003W   | MRPL32    | Stall    |  -0.869|
| YBR201W   | DER1      | No stall |  -0.870|
| YDR343C   | HXT6      | No stall |  -0.874|
| YPL208W   | RKM1      | No stall |  -0.880|
| YGR108W   | CLB1      | No stall |  -0.883|
| YOL017W   | ESC8      | Stall    |  -0.886|
| YPR047W   | MSF1      | No stall |  -0.893|
| YOR036W   | PEP12     | No stall |  -0.894|
| YGL005C   | COG7      | No stall |  -0.904|
| YGR174C   | CBP4      | No stall |  -0.905|
| YJR057W   | CDC8      | No stall |  -0.908|
| YML129C   | COX14     | No stall |  -0.908|
| YDL142C   | CRD1      | No stall |  -0.910|
| YFR042W   | KEG1      | No stall |  -0.910|
| YOR215C   | AIM41     | No stall |  -0.911|
| YKL137W   | CMC1      | No stall |  -0.914|
| YDR050C   | TPI1      | No stall |  -0.917|
| YNL220W   | ADE12     | No stall |  -0.919|
| YMR228W   | MTF1      | No stall |  -0.919|
| YDR116C   | MRPL1     | No stall |  -0.921|
| YHR043C   | DOG2      | No stall |  -0.925|
| YMR193W   | MRPL24    | Stall    |  -0.928|
| YJR080C   | AIM24     | No stall |  -0.929|
| YDR158W   | HOM2      | No stall |  -0.930|
| YPR188C   | MLC2      | No stall |  -0.935|
| YKL216W   | URA1      | No stall |  -0.935|
| YER042W   | MXR1      | No stall |  -0.935|
| YGL031C   | RPL24A    | Stall    |  -0.936|
| YER069W   | ARG5,6    | No stall |  -0.937|
| YOR065W   | CYT1      | No stall |  -0.938|
| YKR083C   | DAD2      | No stall |  -0.939|
| YGR105W   | VMA21     | No stall |  -0.944|
| YNL090W   | RHO2      | No stall |  -0.952|
| YER048W-A | ISD11     | No stall |  -0.955|
| YOL077W-A | ATP19     | No stall |  -0.964|
| YLR061W   | RPL22A    | No stall |  -0.965|
| YDR316W   | OMS1      | Stall    |  -0.973|
| YOL121C   | RPS19A    | No stall |  -0.973|
| YPR100W   | MRPL51    | No stall |  -0.977|
| YMR158W   | MRPS8     | No stall |  -0.978|
| YBL059C-A | CMC2      | Stall    |  -0.979|
| YHR018C   | ARG4      | No stall |  -0.981|
| YIL121W   | QDR2      | No stall |  -0.982|
| YNR001C   | CIT1      | No stall |  -0.982|
| YDR276C   | PMP3      | No stall |  -0.986|
| YDR139C   | RUB1      | No stall |  -0.986|
| YML041C   | VPS71     | No stall |  -0.988|
| YKL006C-A | SFT1      | No stall |  -1.003|
| YKR039W   | GAP1      | No stall |  -1.004|
| Q0250     | COX2      | No stall |  -1.010|
| YJL066C   | MPM1      | No stall |  -1.010|
| YPL132W   | COX11     | Stall    |  -1.011|
| YBR262C   | MIC12     | No stall |  -1.011|
| YPL057C   | SUR1      | No stall |  -1.017|
| YPL087W   | YDC1      | No stall |  -1.017|
| YDR277C   | MTH1      | No stall |  -1.017|
| YMR188C   | MRPS17    | No stall |  -1.027|
| YOR375C   | GDH1      | No stall |  -1.029|
| YPL118W   | MRP51     | No stall |  -1.033|
| YER094C   | PUP3      | No stall |  -1.037|
| YDR225W   | HTA1      | No stall |  -1.041|
| YNL259C   | ATX1      | No stall |  -1.045|
| YER091C   | MET6      | No stall |  -1.050|
| YGR230W   | BNS1      | No stall |  -1.051|
| YGR215W   | RSM27     | Stall    |  -1.052|
| YFR049W   | YMR31     | No stall |  -1.057|
| YNL135C   | FPR1      | No stall |  -1.058|
| YDR033W   | MRH1      | Stall    |  -1.060|
| YJR118C   | ILM1      | No stall |  -1.064|
| YLR287C-A | RPS30A    | Stall    |  -1.080|
| YBL040C   | ERD2      | No stall |  -1.108|
| YDR377W   | ATP17     | No stall |  -1.114|
| YOR369C   | RPS12     | No stall |  -1.138|
| YJR104C   | SOD1      | No stall |  -1.145|
| YPL220W   | RPL1A     | Stall    |  -1.147|
| YDR512C   | EMI1      | No stall |  -1.153|
| YGR008C   | STF2      | No stall |  -1.156|
| YKL001C   | MET14     | No stall |  -1.162|
| YNL037C   | IDH1      | No stall |  -1.169|
| YOL071W   | SDH5      | No stall |  -1.173|
| YLR333C   | RPS25B    | Stall    |  -1.204|
| YIL139C   | REV7      | No stall |  -1.213|
| YOL126C   | MDH2      | No stall |  -1.226|
| YKL003C   | MRP17     | No stall |  -1.228|
| YPL124W   | SPC29     | No stall |  -1.230|
| YBR120C   | CBP6      | No stall |  -1.243|
| YPL013C   | MRPS16    | No stall |  -1.251|
| YMR256C   | COX7      | No stall |  -1.253|
| YKL096W   | CWP1      | No stall |  -1.263|
| Q0275     | COX3      | No stall |  -1.272|
| YNL185C   | MRPL19    | No stall |  -1.303|
| YLR204W   | QRI5      | Stall    |  -1.307|
| YLR325C   | RPL38     | No stall |  -1.311|
| YPR133W-A | TOM5      | No stall |  -1.369|
| YCL025C   | AGP1      | No stall |  -1.406|
| YJR101W   | RSM26     | No stall |  -1.408|
| YKL109W   | HAP4      | Stall    |  -1.420|
| YKL067W   | YNK1      | No stall |  -1.489|
| YML009C   | MRPL39    | Stall    |  -1.491|
| YNL015W   | PBI2      | No stall |  -1.508|
| YHR053C   | CUP1-1    | No stall |  -1.533|
| YHR055C   | CUP1-1    | No stall |  -1.549|
| YDL137W   | ARF2      | No stall |  -1.573|
| YGL227W   | VID30     | No stall |  -1.622|
| YGL255W   | ZRT1      | No stall |  -1.713|
