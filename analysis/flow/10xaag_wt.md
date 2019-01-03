Fluorescence of PGK1 WT reporters with 10xAAG / 10xAGA inserts and varying 5'UTR mutations
================
rasi
02 January, 2019

-   [Import libraries and analysis specific parameters](#import-libraries-and-analysis-specific-parameters)
-   [Read data](#read-data)
-   [Read annotations](#read-annotations)
-   [Rename and calculate average values of fluorescence channels in each well](#rename-and-calculate-average-values-of-fluorescence-channels-in-each-well)
-   [Calculate mean and standard error over replicates](#calculate-mean-and-standard-error-over-replicates)
-   [Plot and tabulate background subtracted and normalized YFP/RFP ratio as a function of initiation codon](#plot-and-tabulate-background-subtracted-and-normalized-yfprfp-ratio-as-a-function-of-initiation-codon)

Import libraries and analysis specific parameters
=================================================

``` r
# standard analysis and plotting functions, includes dplyr, ggplot2 
library(tidyverse)
# loads lab default ggplot2 theme and provides color-blind friendly palette
library(rasilabRtemplates)
# standard error
library(plotrix)

# initiation sites are arranged in this order
initiationmutation_order <- seq(1,8)
names(initiationmutation_order) <- toupper(c( 'ctgc', 'ccgc', 
                              'acgc', 'ccga', 'ccac', 'ccaa', 'caaa', 'aaaa'))

# this folder contains the data and annotations
fcs_file_folder <- "../../data/flow/10xaag_wt/"
```

Read data
=========

``` r
flowdata  <- read_tsv(paste0(fcs_file_folder, '/data.tsv')) %>% 
  print()
```

    ## # A tibble: 660,000 x 7
    ##    plate well   FSC.A  SSC.A FITC.A PE.Texas.Red.A  Time
    ##    <int> <chr>  <int>  <int>  <int>          <int> <dbl>
    ##  1     1 B2     41718  37703     57             41  3.04
    ##  2     1 B2     20876  14209      7            105  3.04
    ##  3     1 B2     34889  22390     56            194  3.04
    ##  4     1 B2     44640  49289    168             15  3.05
    ##  5     1 B2     30240  35159     79             23  3.06
    ##  6     1 B2    130783 109650    333            131  3.08
    ##  7     1 B2     51243  49906     63             19  3.12
    ##  8     1 B2     41055  45669     79            -10  3.14
    ##  9     1 B2     40873  35228     82            -31  3.16
    ## 10     1 B2     36060  27464    131             87  3.16
    ## # ... with 659,990 more rows

Read annotations
================

``` r
annotations  <- read_tsv(paste0(fcs_file_folder, '/annotations.tsv')) %>% 
  print()
```

    ## # A tibble: 132 x 7
    ##    plate well  strain  replicate initiationmutation codonmutation gene   
    ##    <int> <chr> <chr>       <int> <chr>              <chr>         <chr>  
    ##  1     1 B2    by4741          1 CAAA               <NA>          <NA>   
    ##  2     1 B3    schp15          1 CAAA               <NA>          <NA>   
    ##  3     1 B4    schp19          1 CAAA               cgg           maxhis3
    ##  4     1 B5    schp20          1 CAAA               aga           maxhis3
    ##  5     1 B7    schp674         1 CAAA               aag           pgk1   
    ##  6     1 B8    schp675         1 CCGC               aag           pgk1   
    ##  7     1 B9    schp676         1 CCAA               aag           pgk1   
    ##  8     1 B10   schp677         1 CCAC               aag           pgk1   
    ##  9     1 B11   schp678         1 CCGA               aag           pgk1   
    ## 10     1 C2    schp679         1 CTGC               aag           pgk1   
    ## # ... with 122 more rows

Rename and calculate average values of fluorescence channels in each well
=========================================================================

``` r
by_file <- flowdata  %>% 
  # group by each plate and well
  group_by(plate, well) %>% 
  select(FITC.A, PE.Texas.Red.A) %>% 
  # calculate mean
  summarise_all(mean) %>% 
  # rename
  rename('yfp' = FITC.A, 'rfp' = PE.Texas.Red.A) %>% 
  # join annotations
  left_join(annotations, by = c('plate', 'well')) %>% 
  print()
```

    ## # A tibble: 66 x 9
    ## # Groups:   plate [?]
    ##    plate well      yfp     rfp strain  replicate initiationmutation
    ##    <int> <chr>   <dbl>   <dbl> <chr>       <int> <chr>             
    ##  1     1 B10    4379.  20725.  schp677         1 CCAC              
    ##  2     1 B11    4806.  21500.  schp678         1 CCGA              
    ##  3     1 B2       66.4    28.5 by4741          1 CAAA              
    ##  4     1 B3      267.  21366.  schp15          1 CAAA              
    ##  5     1 B4     2327.  18761.  schp19          1 CAAA              
    ##  6     1 B5    20649.  18753.  schp20          1 CAAA              
    ##  7     1 B7     1657.  20194.  schp674         1 CAAA              
    ##  8     1 B8     4136.  21386.  schp675         1 CCGC              
    ##  9     1 B9     3119.  20731.  schp676         1 CCAA              
    ## 10     1 C10    5023.  21059.  schp687         1 CCGA              
    ## # ... with 56 more rows, and 2 more variables: codonmutation <chr>,
    ## #   gene <chr>

Calculate mean and standard error over replicates
=================================================

``` r
avg_data  <- by_file %>% 
  # anti_join(bad_wells) %>% 
  # strain is used to get replicates
  group_by(strain) %>% 
  # calculate mean and std.err
  mutate(mean_yfp = mean(yfp), 
         mean_rfp = mean(rfp)) %>% 
  ungroup() %>% 
  mutate(yfp = yfp - mean_yfp[strain == "schp15" & replicate == 1], 
         rfp = rfp - mean_rfp[strain == "by4741" & replicate == 1]) %>% 
  mutate(yfp_rfp_ratio = yfp / rfp) %>% 
  # calculate mean and standard error
  group_by(strain) %>% 
  mutate(mean_yfp = mean(yfp), 
         mean_rfp = mean(rfp), 
         mean_ratio = mean(yfp_rfp_ratio), 
         se_yfp = std.error(yfp), 
         se_rfp = std.error(rfp),
         se_ratio = std.error(yfp_rfp_ratio),
         n = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  print()
```

    ## # A tibble: 22 x 17
    ##    plate well        yfp       rfp strain  replicate initiationmutation
    ##    <int> <chr>     <dbl>     <dbl> <chr>       <int> <chr>             
    ##  1     1 B2      -198.       -4.58 by4741          1 CAAA              
    ##  2     1 B3         2.93  21333.   schp15          1 CAAA              
    ##  3     1 B4      2063.    18728.   schp19          1 CAAA              
    ##  4     1 B5     20385.    18720.   schp20          1 CAAA              
    ##  5     1 B7      1393.    20161.   schp674         1 CAAA              
    ##  6     1 B8      3872.    21353.   schp675         1 CCGC              
    ##  7     1 B9      2855.    20698.   schp676         1 CCAA              
    ##  8     1 B10     4115.    20691.   schp677         1 CCAC              
    ##  9     1 B11     4542.    21467.   schp678         1 CCGA              
    ## 10     1 C2      2597.    20870.   schp679         1 CTGC              
    ## # ... with 12 more rows, and 10 more variables: codonmutation <chr>,
    ## #   gene <chr>, mean_yfp <dbl>, mean_rfp <dbl>, yfp_rfp_ratio <dbl>,
    ## #   mean_ratio <dbl>, se_yfp <dbl>, se_rfp <dbl>, se_ratio <dbl>, n <int>

``` r
normalization <- avg_data %>% 
  filter(strain == "schp19")
```

Plot and tabulate background subtracted and normalized YFP/RFP ratio as a function of initiation codon
======================================================================================================

``` r
plot_data <- avg_data %>% 
  mutate(mean_ratio = mean_ratio / normalization[[1, "mean_ratio"]]) %>% 
  mutate(se_ratio = se_ratio / normalization[[1, "mean_ratio"]]) %>% 
  filter(initiationmutation != "CTG") %>%
  # arrange initiationmutation in this order
  mutate(initiationmutation = fct_reorder(
      initiationmutation,
      initiationmutation_order[initiationmutation])) %>%
  filter(gene == "pgk1") %>% 
  mutate(codonmutation = paste0("10×", toupper(codonmutation)))

plot_data %>% 
  ggplot(aes(x = initiationmutation, y = mean_ratio, 
             ymin = mean_ratio - se_ratio, ymax = mean_ratio + se_ratio,
             group = codonmutation)) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_errorbar(width = 0.75) +
  facet_wrap(~fct_rev(codonmutation), ncol = 1, scales = "free") + 
  labs(y = 'fluorescence (a.u.)',
       x = '-4 to -1 nt from ATG') +
  theme(legend.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=4))
```

![](10xaag_wt_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
ggsave('figures/initiation_pgk1_aag_flow.pdf')

plot_data %>% 
  arrange(codonmutation, initiationmutation) %>% 
  select(codonmutation, initiationmutation, mean_ratio, se_ratio, n) %>% 
  mutate_if(is.numeric, funs(round(., 3))) %>% 
  knitr::kable()
```

| codonmutation | initiationmutation |  mean\_ratio|  se\_ratio|    n|
|:--------------|:-------------------|------------:|----------:|----:|
| 10×AAG        | CTGC               |        1.109|      0.032|    3|
| 10×AAG        | CCGC               |        1.616|      0.027|    3|
| 10×AAG        | ACGC               |        1.681|      0.033|    3|
| 10×AAG        | CCGA               |        1.920|      0.006|    3|
| 10×AAG        | CCAC               |        1.771|      0.022|    3|
| 10×AAG        | CCAA               |        1.231|      0.026|    3|
| 10×AAG        | CAAA               |        0.632|      0.037|    3|
| 10×AAG        | AAAA               |        0.519|      0.018|    3|
| 10×AGA        | CTGC               |        0.761|      0.107|    3|
| 10×AGA        | CCGC               |        1.446|      0.028|    3|
| 10×AGA        | ACGC               |        1.572|      0.027|    3|
| 10×AGA        | CCGA               |        2.050|      0.007|    3|
| 10×AGA        | CCAC               |        2.144|      0.022|    3|
| 10×AGA        | CCAA               |        2.881|      0.072|    3|
| 10×AGA        | CAAA               |        3.564|      0.158|    3|
| 10×AGA        | AAAA               |        3.638|      0.029|    3|
