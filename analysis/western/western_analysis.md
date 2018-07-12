Analysis of Western blot quantification from PGK1 reporters in different genetic backgrounds
================
rasi
03 January, 2019

-   [Import libraries and analysis-specific parameters](#import-libraries-and-analysis-specific-parameters)
-   [Read data and annotations](#read-data-and-annotations)
-   [Group by each blot, subtract background, normalize FLAG by H3 for each lane, and then by max ratio](#group-by-each-blot-subtract-background-normalize-flag-by-h3-for-each-lane-and-then-by-max-ratio)

Import libraries and analysis-specific parameters
=================================================

``` r
# standard analysis and plotting functions, includes dplyr, ggplot2 
library(tidyverse)
# loads lab default ggplot2 theme and provides color-blind friendly palette
library(rasilabRtemplates)

# initiation sites are arranged in this order
initiationmutation_order <- seq(1,3)
names(initiationmutation_order) <- toupper(c( 'ctgc', 'ccac', 'caaa'))

# this folder contains the data and annotations
data_folder <- "../../data/western/"
```

Read data and annotations
=========================

``` r
data  <- read_tsv(paste0(data_folder, '/quantification.tsv')) %>% 
  print()
```

    ## # A tibble: 112 x 10
    ##     lane  area  mean   min   max codon initiation  blot strain antibody
    ##    <int> <int> <dbl> <int> <int> <chr> <chr>      <int> <chr>  <chr>   
    ##  1     1  2640  179.   123   255 aga   ctgc           1 wt     h3      
    ##  2     2  2640  175.   123   255 aga   ccac           1 wt     h3      
    ##  3     3  2640  177.   124   255 aga   caaa           1 wt     h3      
    ##  4     4  2640  176.   123   255 aag   ctgc           1 wt     h3      
    ##  5     5  2640  174.   124   255 aag   ccac           1 wt     h3      
    ##  6     6  2640  171.   123   255 aag   caaa           1 wt     h3      
    ##  7     7  2640  124.   122   129 <NA>  <NA>           1 wt     h3      
    ##  8     1  4000  134.   100   244 aga   ctgc           1 wt     flag    
    ##  9     2  4000  173.    99   255 aga   ccac           1 wt     flag    
    ## 10     3  4000  208.    99   255 aga   caaa           1 wt     flag    
    ## # ... with 102 more rows

Group by each blot, subtract background, normalize FLAG by H3 for each lane, and then by max ratio
==================================================================================================

``` r
data %>% 
  group_by(blot, antibody) %>% 
  # subtract background
  mutate(mean = mean - mean[is.na(codon) & is.na(initiation)]) %>% 
  ungroup() %>% 
  # get rid of unwanted columns
  select(-area, -min, -max) %>% 
  # get rid of background lanes
  filter(!is.na(codon) & !is.na(initiation)) %>% 
  # get FLAG and H3 quantification side by side for each lane
  spread(antibody, mean) %>% 
  # normalize FLAG by H3
  mutate(ratio = flag / h3) %>% 
  select(-flag, -h3) %>% 
  # normalize by maximum ratio within each blot
  group_by(blot) %>% 
  mutate(ratio = ratio / max(ratio)) %>%
  ungroup() %>% 
  # multiply by 10 and convert to integer for ease of display
  mutate(ratio = as.integer(ratio * 10)) %>%
  arrange(blot, lane) %>% 
  select(blot, lane, strain, codon, initiation, ratio) %>%
  knitr::kable()
```

|  blot|  lane| strain | codon | initiation |  ratio|
|-----:|-----:|:-------|:------|:-----------|------:|
|     1|     1| wt     | aga   | ctgc       |      3|
|     1|     2| wt     | aga   | ccac       |      7|
|     1|     3| wt     | aga   | caaa       |     10|
|     1|     4| wt     | aag   | ctgc       |      2|
|     1|     5| wt     | aag   | ccac       |      6|
|     1|     6| wt     | aag   | caaa       |      4|
|     2|     1| wt     | cca   | ctgc       |      7|
|     2|     2| wt     | cca   | ccac       |      8|
|     2|     3| wt     | cca   | caaa       |     10|
|     2|     4| wt     | ccg   | ctgc       |      6|
|     2|     5| wt     | ccg   | ccac       |      7|
|     2|     6| wt     | ccg   | caaa       |      4|
|     3|     1| wt     | aga   | ctgc       |      7|
|     3|     2| wt     | aga   | ccac       |      8|
|     3|     3| wt     | aga   | caaa       |     10|
|     3|     4| wt     | cgg   | ctgc       |      7|
|     3|     5| wt     | cgg   | ccac       |      7|
|     3|     6| wt     | cgg   | caaa       |      6|
|     4|     1| wt     | ccg   | ctgc       |      5|
|     4|     2| wt     | ccg   | ccac       |      6|
|     4|     3| wt     | ccg   | caaa       |      5|
|     4|     4| wt     | cca   | ctgc       |      5|
|     4|     5| wt     | cca   | ccac       |      7|
|     4|     6| wt     | cca   | caaa       |     10|
|     5|     1| ltn1   | ccg   | ctgc       |      7|
|     5|     2| ltn1   | ccg   | ccac       |      9|
|     5|     3| ltn1   | ccg   | caaa       |      7|
|     5|     4| ltn1   | cca   | ctgc       |      6|
|     5|     5| ltn1   | cca   | ccac       |      7|
|     5|     6| ltn1   | cca   | caaa       |     10|
|     6|     1| dom34  | ccg   | ctgc       |      5|
|     6|     2| dom34  | ccg   | ccac       |      7|
|     6|     3| dom34  | ccg   | caaa       |      5|
|     6|     4| dom34  | cca   | ctgc       |      5|
|     6|     5| dom34  | cca   | ccac       |      8|
|     6|     6| dom34  | cca   | caaa       |     10|
|     7|     1| hel2   | ccg   | ctgc       |      4|
|     7|     2| hel2   | ccg   | ccac       |      6|
|     7|     3| hel2   | ccg   | caaa       |      6|
|     7|     4| hel2   | cca   | ctgc       |      5|
|     7|     5| hel2   | cca   | ccac       |      7|
|     7|     6| hel2   | cca   | caaa       |     10|
|     8|     1| asc1   | ccg   | ctgc       |      4|
|     8|     2| asc1   | ccg   | ccac       |      7|
|     8|     3| asc1   | ccg   | caaa       |     10|
|     8|     4| asc1   | cca   | ctgc       |      4|
|     8|     5| asc1   | cca   | ccac       |      4|
|     8|     6| asc1   | cca   | caaa       |      6|
