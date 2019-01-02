Analyze simulation results
================
rasi
01 January, 2019

-   [Load libraries](#load-libraries)
-   [Read protein count data](#read-protein-count-data)
-   [Read mRNA lifetime data](#read-mrna-lifetime-data)
-   [Read simulation parameters](#read-simulation-parameters)
-   [Combine all data into a single table](#combine-all-data-into-a-single-table)
-   [Proteins per mRNA as a function of initiation rate](#proteins-per-mrna-as-a-function-of-initiation-rate)
-   [mRNA lifetime as a function of initiation rate](#mrna-lifetime-as-a-function-of-initiation-rate)

Load libraries
==============

``` r
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=999)
```

Read protein count data
=======================

``` r
psr_data <- read_tsv("tables/psr_stats.tsv") %>% 
  print()
```

    ## # A tibble: 9 x 6
    ##   sim_id mean_p_per_m sd_p_per_m total_p total_time     psr
    ##    <int>        <int>      <int>   <int>      <int>   <dbl>
    ## 1      0            3          2    3416     999483 0.00342
    ## 2      1            7          4    6909     999516 0.00691
    ## 3      2           14          7   13347     999970 0.0133 
    ## 4      3           26         13   27232     999067 0.0273 
    ## 5      4           54         25   55366     999954 0.0554 
    ## 6      5          109         51  104055     999998 0.104  
    ## 7      6          208         86  198587     999999 0.199  
    ## 8      7          388        133  384541     999999 0.385  
    ## 9      8          672        184  235298     359257 0.655

Read mRNA lifetime data
=======================

``` r
lifetime_data <- read_tsv("tables/mrna_lifetime_stats.tsv") %>% 
  mutate(se_lifetime = sd_lifetime / sqrt(n_mrna)) %>% 
  print()
```

    ## # A tibble: 9 x 5
    ##   sim_id mean_lifetime sd_lifetime n_mrna se_lifetime
    ##    <int>         <int>       <int>  <int>       <dbl>
    ## 1      0          2108         271    965        8.72
    ## 2      1          2112         277    969        8.90
    ## 3      2          2110         277    938        9.04
    ## 4      3          2089         284   1032        8.84
    ## 5      4          2108         265   1006        8.36
    ## 6      5          2085         280    952        9.07
    ## 7      6          2086         276    951        8.95
    ## 8      7          2106         288    987        9.17
    ## 9      8          2095         273    347       14.7

Read simulation parameters
==========================

``` r
annotations <- list.files("output/", pattern = "params.tsv.gz$", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(sim_id = str_extract(file, "(?<=tasep_)[[:digit:]]+")) %>% 
  mutate(data = map(file, read_tsv)) %>% 
  select(-file, -sno) %>% 
  unnest() %>% 
  type_convert() %>% 
  filter(parameter == 'k_init') %>% 
  spread(parameter, value) %>% 
  print()
```

    ## # A tibble: 9 x 2
    ##   sim_id  k_init
    ##    <int>   <dbl>
    ## 1      0 0.00391
    ## 2      1 0.00781
    ## 3      2 0.0156 
    ## 4      3 0.0312 
    ## 5      4 0.0625 
    ## 6      5 0.125  
    ## 7      6 0.250  
    ## 8      7 0.500  
    ## 9      8 1.00

Combine all data into a single table
====================================

``` r
data <- annotations %>% 
  left_join(psr_data, by = "sim_id") %>% 
  left_join(lifetime_data, by = "sim_id") %>% 
  print()
```

    ## # A tibble: 9 x 11
    ##   sim_id  k_init mean_p_per_m sd_p_per_m total_p total_time     psr
    ##    <int>   <dbl>        <int>      <int>   <int>      <int>   <dbl>
    ## 1      0 0.00391            3          2    3416     999483 0.00342
    ## 2      1 0.00781            7          4    6909     999516 0.00691
    ## 3      2 0.0156            14          7   13347     999970 0.0133 
    ## 4      3 0.0312            26         13   27232     999067 0.0273 
    ## 5      4 0.0625            54         25   55366     999954 0.0554 
    ## 6      5 0.125            109         51  104055     999998 0.104  
    ## 7      6 0.250            208         86  198587     999999 0.199  
    ## 8      7 0.500            388        133  384541     999999 0.385  
    ## 9      8 1.00             672        184  235298     359257 0.655  
    ## # ... with 4 more variables: mean_lifetime <int>, sd_lifetime <int>,
    ## #   n_mrna <int>, se_lifetime <dbl>

Proteins per mRNA as a function of initiation rate
==================================================

``` r
plot_data <- data

plot_data %>%
  ggplot(aes(x = k_init, y = mean_p_per_m)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_p_per_m - sd_p_per_m, 
                    ymax = mean_p_per_m + sd_p_per_m, width = 0.5)) +
  scale_x_continuous(trans = "log2",
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     breaks = 2^(seq(-8,0,2))) +
  scale_y_continuous(trans = "log10",
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     breaks = 10^(seq(1,3))) +
  # scale_y_continuous(limits = c(0, NA)) +
  labs(x = "initiation rate (s-1)", y = "proteins per mRNA") +
  theme(legend.position = "top")
```

![](/fh/fast/subramaniam_a/user/rasi/git/ribosome_collisions_yeast/modeling/simulation_runs/tasep_with_canonical_decay/analyze_results_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
ggsave("figures/protein_per_mrna_vs_initiation_rate.pdf")
```

mRNA lifetime as a function of initiation rate
==============================================

``` r
plot_data <- data

plot_data %>%
  ggplot(aes(x = k_init, y = mean_lifetime)) +
  geom_point(size = 2, shape = 17) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_lifetime - sd_lifetime, 
                    ymax = mean_lifetime + sd_lifetime, width = 0.5)) +
  scale_x_continuous(trans = "log2",
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     breaks = 2^(seq(-8,0,2))) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "initiation rate (s-1)", y = "mean mRNA lifetime (s)") +
  theme(legend.position = "top")
```

![](/fh/fast/subramaniam_a/user/rasi/git/ribosome_collisions_yeast/modeling/simulation_runs/tasep_with_canonical_decay/analyze_results_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
ggsave("figures/mrna_lifetime_vs_initiation_rate.pdf")
```
