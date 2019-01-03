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

    ## # A tibble: 150,000 x 7
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
    ## # ... with 149,990 more rows

Read annotations
================

``` r
annotations  <- read_tsv(paste0(fcs_file_folder, '/annotations.tsv')) %>% 
  print()
```

    ## # A tibble: 30 x 8
    ##    plate well  strain replicate initiationmutation codonmutation gene   
    ##    <int> <chr> <chr>      <int> <chr>              <chr>         <chr>  
    ##  1     1 B2    by4741         1 CAAA               <NA>          <NA>   
    ##  2     1 B3    schp15         1 CAAA               <NA>          <NA>   
    ##  3     1 B4    schp19         1 CAAA               cgg           maxhis3
    ##  4     1 B5    schp20         1 CAAA               aga           maxhis3
    ##  5     1 B6    schp76         1 CAAA               aga           pgk1   
    ##  6     1 E2    by4741         2 CAAA               <NA>          <NA>   
    ##  7     1 E3    schp15         2 CAAA               <NA>          <NA>   
    ##  8     1 E4    schp19         2 CAAA               cgg           maxhis3
    ##  9     1 E5    schp20         2 CAAA               aga           maxhis3
    ## 10     1 E6    schp76         2 CAAA               aga           pgk1   
    ## # ... with 20 more rows, and 1 more variable: knockout <chr>

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

    ## # A tibble: 15 x 10
    ## # Groups:   plate [?]
    ##    plate well      yfp     rfp strain replicate initiationmutation
    ##    <int> <chr>   <dbl>   <dbl> <chr>      <int> <chr>             
    ##  1     1 B2       66.4    28.5 by4741         1 CAAA              
    ##  2     1 B3      267.  21366.  schp15         1 CAAA              
    ##  3     1 B4     2327.  18761.  schp19         1 CAAA              
    ##  4     1 B5    20649.  18753.  schp20         1 CAAA              
    ##  5     1 B6    13230.  23298.  schp76         1 CAAA              
    ##  6     2 E2       64.3    36.9 by4741         4 CAAA              
    ##  7     2 E3      275.  21173.  schp15         4 CAAA              
    ##  8     2 E4     2250.  18882.  schp19         4 CAAA              
    ##  9     2 E5    19844.  19462.  schp20         4 CAAA              
    ## 10     2 E6    12196.  22131.  schp76         4 CAAA              
    ## 11     3 B2       64.4    33.9 by4741         5 CAAA              
    ## 12     3 B3      250.  18786.  schp15         5 CAAA              
    ## 13     3 B4     2224.  17281.  schp19         5 CAAA              
    ## 14     3 B5    18249.  17130.  schp20         5 CAAA              
    ## 15     3 B6    12354.  21320.  schp76         5 CAAA              
    ## # ... with 3 more variables: codonmutation <chr>, gene <chr>,
    ## #   knockout <chr>

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

    ## # A tibble: 5 x 18
    ##   plate well        yfp       rfp strain replicate initiationmutation
    ##   <int> <chr>     <dbl>     <dbl> <chr>      <int> <chr>             
    ## 1     1 B2      -198.       -4.58 by4741         1 CAAA              
    ## 2     1 B3         2.93  21333.   schp15         1 CAAA              
    ## 3     1 B4      2063.    18728.   schp19         1 CAAA              
    ## 4     1 B5     20385.    18720.   schp20         1 CAAA              
    ## 5     1 B6     12966.    23265.   schp76         1 CAAA              
    ## # ... with 11 more variables: codonmutation <chr>, gene <chr>,
    ## #   knockout <chr>, mean_yfp <dbl>, mean_rfp <dbl>, yfp_rfp_ratio <dbl>,
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
| 10×AGA        | CAAA               |         5.06|      0.074|    3|
