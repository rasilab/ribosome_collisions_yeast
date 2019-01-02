Analyze simulation results
================
rasi
01 January, 2019

-   [Load libraries](#load-libraries)
-   [Read protein count data](#read-protein-count-data)
-   [Read simulation parameters](#read-simulation-parameters)
-   [Combine all data into a single table](#combine-all-data-into-a-single-table)
-   [Protein synthesis rate as a function of abortive termination rate](#protein-synthesis-rate-as-a-function-of-abortive-termination-rate)

Load libraries
--------------

``` r
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=999)

preterm_model_names <- c(
  "hit5" = "CSAT",
  "simple" = "SAT",
  "trafficjam" = "TJ"
)
```

Read protein count data
=======================

``` r
psr_data <- read_tsv("tables/psr_stats.tsv") %>% 
  print()
```

    ## # A tibble: 60 x 6
    ##    sim_id mean_p_per_m sd_p_per_m total_p total_time       psr
    ##     <int>        <int> <chr>        <int>      <int>     <dbl>
    ##  1      0        59977 <NA>         59977     657244 0.0913   
    ##  2      1        70962 <NA>         70962     786155 0.0903   
    ##  3     10        14479 <NA>         14479     999997 0.0145   
    ##  4     11         5955 <NA>          5955     999999 0.00596  
    ##  5     12         1745 <NA>          1745     999998 0.00174  
    ##  6     13          312 <NA>           312     999998 0.000312 
    ##  7     14           38 <NA>            38     999993 0.0000380
    ##  8     15        71583 <NA>         71583     999999 0.0716   
    ##  9     16        67402 <NA>         67402     999995 0.0674   
    ## 10     17        62676 <NA>         62676     999999 0.0627   
    ## # ... with 50 more rows

Read simulation parameters
==========================

``` r
sim_params <- read_tsv("sim.params.tsv") %>% 
  rename(sim_id = X1) %>% 
  mutate(k_elong_stall = str_split(k_elong_stall, ",")) %>%
  mutate(k_elong_stall = map(k_elong_stall, as.numeric)) %>%
  mutate(k_elong_stall = map(k_elong_stall, function(x) unique(x))) %>%
  unnest() %>%
  mutate(x_stall = stringr::str_split(x_stall, ',')) %>%
  mutate(k_stall = k_elong_stall / as.numeric(n_stall)) %>%
  select(sim_id, preterm_intact_rate, preterm_intact_model, k_stall) %>% 
  print()
```

    ## # A tibble: 60 x 4
    ##    sim_id preterm_intact_rate preterm_intact_model k_stall
    ##     <int>               <dbl> <chr>                  <dbl>
    ##  1      0            0.000884 simple                 0.100
    ##  2      1            0.00125  simple                 0.100
    ##  3      2            0.00177  simple                 0.100
    ##  4      3            0.00250  simple                 0.100
    ##  5      4            0.00354  simple                 0.100
    ##  6      5            0.00500  simple                 0.100
    ##  7      6            0.00707  simple                 0.100
    ##  8      7            0.0100   simple                 0.100
    ##  9      8            0.0141   simple                 0.100
    ## 10      9            0.0200   simple                 0.100
    ## # ... with 50 more rows

``` r
annotations <- list.files("output/", pattern = "params.tsv.gz$", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(sim_id = str_extract(file, "(?<=tasep_)[[:digit:]]+")) %>% 
  mutate(data = map(file, read_tsv)) %>% 
  select(-file, -sno) %>% 
  unnest() %>% 
  type_convert() %>%
  # retain only parameters that are varied, the others are for checking
  group_by(parameter) %>%
  mutate(vary = if_else(length(unique(value)) > 1, T, F)) %>%
  ungroup() %>%
  filter(vary == T) %>%
  select(-vary) %>%
  spread(parameter, value) %>%
  left_join(sim_params, by = "sim_id") %>% 
  print()
```

    ## # A tibble: 60 x 13
    ##    sim_id k_elong_stall_1 k_elong_stall_2 k_elong_stall_3 k_elong_stall_4
    ##     <int>           <dbl>           <dbl>           <dbl>           <dbl>
    ##  1      0           0.600           0.600           0.600           0.600
    ##  2      1           0.600           0.600           0.600           0.600
    ##  3      2           0.600           0.600           0.600           0.600
    ##  4      3           0.600           0.600           0.600           0.600
    ##  5      4           0.600           0.600           0.600           0.600
    ##  6      5           0.600           0.600           0.600           0.600
    ##  7      6           0.600           0.600           0.600           0.600
    ##  8      7           0.600           0.600           0.600           0.600
    ##  9      8           0.600           0.600           0.600           0.600
    ## 10      9           0.600           0.600           0.600           0.600
    ## # ... with 50 more rows, and 8 more variables: k_elong_stall_5 <dbl>,
    ## #   k_elong_stall_6 <dbl>, k_preterm_5_hit_intact <dbl>,
    ## #   k_preterm_both_hit_intact <dbl>, k_preterm_no_hit_intact <dbl>,
    ## #   preterm_intact_rate <dbl>, preterm_intact_model <chr>, k_stall <dbl>

Combine all data into a single table
====================================

``` r
data <- annotations %>% 
  left_join(psr_data, by = "sim_id") %>% 
  print()
```

    ## # A tibble: 60 x 18
    ##    sim_id k_elong_stall_1 k_elong_stall_2 k_elong_stall_3 k_elong_stall_4
    ##     <int>           <dbl>           <dbl>           <dbl>           <dbl>
    ##  1      0           0.600           0.600           0.600           0.600
    ##  2      1           0.600           0.600           0.600           0.600
    ##  3      2           0.600           0.600           0.600           0.600
    ##  4      3           0.600           0.600           0.600           0.600
    ##  5      4           0.600           0.600           0.600           0.600
    ##  6      5           0.600           0.600           0.600           0.600
    ##  7      6           0.600           0.600           0.600           0.600
    ##  8      7           0.600           0.600           0.600           0.600
    ##  9      8           0.600           0.600           0.600           0.600
    ## 10      9           0.600           0.600           0.600           0.600
    ## # ... with 50 more rows, and 13 more variables: k_elong_stall_5 <dbl>,
    ## #   k_elong_stall_6 <dbl>, k_preterm_5_hit_intact <dbl>,
    ## #   k_preterm_both_hit_intact <dbl>, k_preterm_no_hit_intact <dbl>,
    ## #   preterm_intact_rate <dbl>, preterm_intact_model <chr>, k_stall <dbl>,
    ## #   mean_p_per_m <int>, sd_p_per_m <chr>, total_p <int>, total_time <int>,
    ## #   psr <dbl>

Protein synthesis rate as a function of abortive termination rate
=================================================================

``` r
plot_data <- data %>% 
  mutate(model = preterm_model_names[preterm_intact_model]) %>% 
  mutate(model_stall = paste0(model, "_", k_stall))

plot_data %>%
  ggplot(aes(x = preterm_intact_rate, y = psr, color = model, shape = model_stall)) +
  facet_wrap(~model, ncol = 2, scales = "free") +
  geom_line() +
  geom_point(size = 1.5) +
  scale_x_continuous(trans = "log2",
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     breaks = scales::trans_breaks("log2", function(x) 2^x, n = 5)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = cbPalette[c(3,2)]) +
  scale_shape_manual(values = c(19, 1, 17, 2)) +
  labs(x = "abortive termination rate (s-1)", 
       y = "protein synthesis rate (s-1)", color = "", shape = "") +
  theme(legend.position = "top")
```

![](/fh/fast/subramaniam_a/user/rasi/git/ribosome_collisions_yeast/modeling/simulation_runs/preterm_rate_vary/analyze_results_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
ggsave("figures/psr_vs_abortive_termination_rate.pdf", width = 3, height = 2.2)
```
