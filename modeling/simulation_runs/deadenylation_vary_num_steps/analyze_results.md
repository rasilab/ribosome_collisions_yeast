Analyze simulation results
================
rasi
01 January, 2019

-   [Load libraries](#load-libraries)
-   [Read mRNA lifetime data](#read-mrna-lifetime-data)
-   [Read simulation parameters](#read-simulation-parameters)
-   [Combine all data into a single table](#combine-all-data-into-a-single-table)
-   [mRNA lifetime as a function of number of deadenylation steps](#mrna-lifetime-as-a-function-of-number-of-deadenylation-steps)

Load libraries
==============

``` r
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=999)
```

Read mRNA lifetime data
=======================

``` r
lifetime_data <- read_tsv("tables/mrna_lifetime_stats.tsv") %>% 
  mutate(se_lifetime = sd_lifetime / sqrt(n_mrna)) %>% 
  print()
```

    ## # A tibble: 7 x 5
    ##   sim_id mean_lifetime sd_lifetime n_mrna se_lifetime
    ##    <int>         <int>       <int>  <int>       <dbl>
    ## 1      0          1981        1969    987       62.7 
    ## 2      1          2196        1491    996       47.2 
    ## 3      2          2096         978    963       31.5 
    ## 4      3          2081         698    997       22.1 
    ## 5      4          2090         521   1035       16.2 
    ## 6      5          2111         363    961       11.7 
    ## 7      6          2079         271    963        8.73

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
  filter(parameter %in% c('l_polya', 'k_deadenylation')) %>% 
  spread(parameter, value) %>% 
  print()
```

    ## # A tibble: 7 x 3
    ##   sim_id k_deadenylation l_polya
    ##    <int>           <dbl>   <dbl>
    ## 1      0        0.000500      1.
    ## 2      1        0.00100       2.
    ## 3      2        0.00200       4.
    ## 4      3        0.00400       8.
    ## 5      4        0.00800      16.
    ## 6      5        0.0160       32.
    ## 7      6        0.0320       64.

Combine all data into a single table
====================================

``` r
data <- annotations %>% 
  left_join(lifetime_data, by = "sim_id") %>% 
  print()
```

    ## # A tibble: 7 x 7
    ##   sim_id k_deadenylation l_polya mean_lifetime sd_lifetime n_mrna
    ##    <int>           <dbl>   <dbl>         <int>       <int>  <int>
    ## 1      0        0.000500      1.          1981        1969    987
    ## 2      1        0.00100       2.          2196        1491    996
    ## 3      2        0.00200       4.          2096         978    963
    ## 4      3        0.00400       8.          2081         698    997
    ## 5      4        0.00800      16.          2090         521   1035
    ## 6      5        0.0160       32.          2111         363    961
    ## 7      6        0.0320       64.          2079         271    963
    ## # ... with 1 more variable: se_lifetime <dbl>

mRNA lifetime as a function of number of deadenylation steps
============================================================

``` r
plot_data <- data

plot_data %>%
  ggplot(aes(x = l_polya, y = mean_lifetime)) +
  geom_point(size = 2, shape = 17) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_lifetime - sd_lifetime, 
                    ymax = mean_lifetime + sd_lifetime, width = 0.5)) +
  scale_x_continuous(trans = "log2",
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     breaks = 2^(seq(0,6,2))) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "number of\ndeadenylation steps", y = "mRNA lifetime (s)") +
  theme(legend.position = "top")
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-6-1.png)
![](analyze_results_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
ggsave("figures/mrna_lifetime_vs_l_polya.pdf")
```
