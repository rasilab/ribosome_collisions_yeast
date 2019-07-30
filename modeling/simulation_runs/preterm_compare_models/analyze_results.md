Analyze simulation results
================
rasi
30 July, 2019

-   [Load libraries](#load-libraries)
-   [Read protein count data](#read-protein-count-data)
-   [Read collision data](#read-collision-data)
    -   [Read simulation parameters](#read-simulation-parameters)
    -   [Combine all data into a single table](#combine-all-data-into-a-single-table)
    -   [How does PSR vary as a function of initiation rate in all models with low stall elongation rate and medium preterm\_intact\_rate?](#how-does-psr-vary-as-a-function-of-initiation-rate-in-all-models-with-low-stall-elongation-rate-and-medium-preterm_intact_rate)
-   [Source data for Fig. 3B](#source-data-for-fig.-3b)
-   [Change in maximal protein synthesis as a function of stall elongation rate in TJ model](#change-in-maximal-protein-synthesis-as-a-function-of-stall-elongation-rate-in-tj-model)
-   [Source data for S1 Fig panel A](#source-data-for-s1-fig-panel-a)
-   [Change in protein synthesis as a function of abort rate in SAT model](#change-in-protein-synthesis-as-a-function-of-abort-rate-in-sat-model)
-   [Source data for S1 Fig panel C](#source-data-for-s1-fig-panel-c)
-   [Change in protein synthesis as a function of stall elongation rate in SAT model](#change-in-protein-synthesis-as-a-function-of-stall-elongation-rate-in-sat-model)
-   [Source data for S1 Fig panel B](#source-data-for-s1-fig-panel-b)
-   [Change in protein synthesis as a function of stall elongation rate in CSAT model](#change-in-protein-synthesis-as-a-function-of-stall-elongation-rate-in-csat-model)
-   [Source data for S1 Fig panel D](#source-data-for-s1-fig-panel-d)

Load libraries
--------------

``` r
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=999)

model_names <- c(
  "hit3" = "CAT",
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

    ## # A tibble: 405 x 6
    ##    sim_id mean_p_per_m sd_p_per_m total_p total_time     psr
    ##     <int>        <int> <chr>        <int>      <int>   <dbl>
    ##  1      0         3851 <NA>          3851     999918 0.00385
    ##  2      1         1275 <NA>          1275     999668 0.00128
    ##  3     10         3245 <NA>          3245     999981 0.00324
    ##  4    100         8886 <NA>          8886     999945 0.00889
    ##  5    101         7613 <NA>          7613     999983 0.00761
    ##  6    102         8844 <NA>          8844     999817 0.00885
    ##  7    103         7412 <NA>          7412     999992 0.00741
    ##  8    104         8673 <NA>          8673     999939 0.00867
    ##  9    105        15427 <NA>         15427     999988 0.0154 
    ## 10    106         7207 <NA>          7207     999956 0.00721
    ## # ... with 395 more rows

Read collision data
===================

``` r
collision_data <- read_tsv("tables/collision_stats.tsv") %>% 
  print()
```

    ## # A tibble: 405 x 6
    ##    sim_id mean_p_per_m sd_p_per_m total_collision total_time collision_freq
    ##     <int>        <int> <chr>                <int>      <int>          <dbl>
    ##  1      0         3754 <NA>                  3754     999918       0.00375 
    ##  2      1         1258 <NA>                  1258     999668       0.00126 
    ##  3     10          965 <NA>                   965     999981       0.000965
    ##  4    100        11022 <NA>                 11022     999945       0.0110  
    ##  5    101        11178 <NA>                 11178     999983       0.0112  
    ##  6    102         9936 <NA>                  9936     999817       0.00994 
    ##  7    103         9298 <NA>                  9298     999992       0.00930 
    ##  8    104         7737 <NA>                  7737     999939       0.00774 
    ##  9    105        15350 <NA>                 15350     999988       0.0154  
    ## 10    106         7637 <NA>                  7637     999956       0.00764 
    ## # ... with 395 more rows

Read simulation parameters
--------------------------

``` r
annotations  <- read_tsv('sim.params.tsv', 
                         col_types = cols(x_stall = col_character(),
                                          k_elong_stall = col_character())) %>%
  mutate(preterm_intact_model = if_else(preterm_intact_rate == 0, 
                                        "trafficjam", 
                                        preterm_intact_model)) %>% 
  rename(sim_id = X1) %>%
  mutate(k_elong_stall = str_split(k_elong_stall, ",")) %>%
  mutate(k_elong_stall = map(k_elong_stall, as.numeric)) %>%
  mutate(k_elong_stall = map(k_elong_stall, function(x) unique(x))) %>%
  unnest() %>%
  mutate(x_stall = stringr::str_split(x_stall, ',')) %>%
  mutate(k_stall = k_elong_stall / as.numeric(n_stall)) %>%
  mutate(n_stall = factor(n_stall)) %>%
  select(sim_id, k_init, k_elong_stall, k_stall, x_stall, n_stall,
         preterm_intact_model, preterm_intact_rate) %>%
  print()
```

    ## # A tibble: 405 x 8
    ##    sim_id  k_init k_elong_stall k_stall x_stall n_stall preterm_intact_…
    ##     <int>   <dbl>         <dbl>   <dbl> <list>  <fct>   <chr>           
    ##  1      0 0.00391          0.12    0.02 <chr [… 6       trafficjam      
    ##  2      1 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  3      2 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  4      3 0.00391          0.12    0.02 <chr [… 6       hit5            
    ##  5      4 0.00391          0.12    0.02 <chr [… 6       hit3            
    ##  6      5 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  7      6 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  8      7 0.00391          0.12    0.02 <chr [… 6       hit5            
    ##  9      8 0.00391          0.12    0.02 <chr [… 6       hit3            
    ## 10      9 0.00391          0.12    0.02 <chr [… 6       hit5            
    ## # ... with 395 more rows, and 1 more variable: preterm_intact_rate <dbl>

Combine all data into a single table
------------------------------------

``` r
data <- annotations %>% 
  left_join(psr_data, by = "sim_id") %>% 
  left_join(collision_data, by = "sim_id") %>% 
  print()
```

    ## # A tibble: 405 x 18
    ##    sim_id  k_init k_elong_stall k_stall x_stall n_stall preterm_intact_…
    ##     <int>   <dbl>         <dbl>   <dbl> <list>  <fct>   <chr>           
    ##  1      0 0.00391          0.12    0.02 <chr [… 6       trafficjam      
    ##  2      1 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  3      2 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  4      3 0.00391          0.12    0.02 <chr [… 6       hit5            
    ##  5      4 0.00391          0.12    0.02 <chr [… 6       hit3            
    ##  6      5 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  7      6 0.00391          0.12    0.02 <chr [… 6       simple          
    ##  8      7 0.00391          0.12    0.02 <chr [… 6       hit5            
    ##  9      8 0.00391          0.12    0.02 <chr [… 6       hit3            
    ## 10      9 0.00391          0.12    0.02 <chr [… 6       hit5            
    ## # ... with 395 more rows, and 11 more variables:
    ## #   preterm_intact_rate <dbl>, mean_p_per_m.x <int>, sd_p_per_m.x <chr>,
    ## #   total_p <int>, total_time.x <int>, psr <dbl>, mean_p_per_m.y <int>,
    ## #   sd_p_per_m.y <chr>, total_collision <int>, total_time.y <int>,
    ## #   collision_freq <dbl>

How does PSR vary as a function of initiation rate in all models with low stall elongation rate and medium preterm\_intact\_rate?
---------------------------------------------------------------------------------------------------------------------------------

``` r
plot_data <- data %>% 
  mutate(k_stall = round(k_stall, 1)) %>%
  filter(k_stall == 0.1) %>%
  filter(((preterm_intact_rate == 0 | preterm_intact_rate == 0.02) 
          & preterm_intact_model == "simple") |
        ((preterm_intact_rate == 0 | preterm_intact_rate == 1) 
                & preterm_intact_model != "simple")) %>% 
  mutate(model = forcats::fct_rev(model_names[preterm_intact_model]))

plot_data %>% 
  ggplot(aes(x = k_init, y = psr, shape = model, fill = model, color = model)) +
  geom_point(size = 1.5) + geom_line(size = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8, 2^0)) +
  scale_y_continuous(breaks = seq(0,0.09, 0.03)) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(25, 24, 21, 23)) +                                                  
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       fill = "", shape = "", color = "") +
  theme(legend.key.height = unit(0.2, "in"))                     
```

![](analyze_results_files/figure-markdown_github/psr_no_quality_control-1.png)

``` r
ggsave('figures/psr_all_models_medium_stall_medium_pretermintact.pdf') 
```

Source data for Fig. 3B
=======================

``` r
plot_data %>% 
  select(model, k_init, psr) %>% 
  knitr::kable()
```

| model |    k\_init|       psr|
|:------|----------:|---------:|
| TJ    |  0.0039062|  0.003810|
| SAT   |  0.0039062|  0.000929|
| CSAT  |  0.0039062|  0.003779|
| CAT   |  0.0039062|  0.003758|
| TJ    |  0.0078125|  0.007653|
| SAT   |  0.0078125|  0.001810|
| CSAT  |  0.0078125|  0.007300|
| CAT   |  0.0078125|  0.007278|
| TJ    |  0.0156250|  0.015427|
| SAT   |  0.0156250|  0.003433|
| CSAT  |  0.0156250|  0.013627|
| CAT   |  0.0156250|  0.013560|
| TJ    |  0.0312500|  0.030367|
| SAT   |  0.0312500|  0.006784|
| CSAT  |  0.0312500|  0.023620|
| CAT   |  0.0312500|  0.024076|
| TJ    |  0.0625000|  0.059270|
| SAT   |  0.0625000|  0.013176|
| CSAT  |  0.0625000|  0.035599|
| CAT   |  0.0625000|  0.038994|
| TJ    |  0.1250000|  0.093576|
| SAT   |  0.1250000|  0.025729|
| CSAT  |  0.1250000|  0.041983|
| CAT   |  0.1250000|  0.055785|
| TJ    |  0.2500000|  0.093712|
| SAT   |  0.2500000|  0.046755|
| CSAT  |  0.2500000|  0.035981|
| CAT   |  0.2500000|  0.069958|
| TJ    |  0.5000000|  0.093668|
| SAT   |  0.5000000|  0.051756|
| CSAT  |  0.5000000|  0.024722|
| CAT   |  0.5000000|  0.078337|
| TJ    |  1.0000000|  0.093551|
| SAT   |  1.0000000|  0.052107|
| CSAT  |  1.0000000|  0.019194|
| CAT   |  1.0000000|  0.081678|

Change in maximal protein synthesis as a function of stall elongation rate in TJ model
======================================================================================

``` r
plot_data <- data %>% 
  filter(preterm_intact_rate == 0)

plot_data %>% 
  ggplot(aes(x = k_init, y = psr, color = as.factor(k_stall), shape = as.factor(k_stall))) +
  geom_point(size = 1.5) + geom_line(size = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8.5, 2^0)) +
  scale_y_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8.5, 2^0)) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(19,17,15)) +
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       shape = "kstall (s-1)", color = "kstall (s-1)") +
  theme(legend.key.height = unit(0.2, "in")) +
  geom_vline(aes(xintercept = k_stall, color = as.factor(k_stall)), 
             show.legend = F, linetype = "dotted")
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
ggsave('figures/psr_tj_model_vary_stall_strength.pdf') 
```

Source data for S1 Fig panel A
==============================

``` r
plot_data %>% 
  select(k_stall, k_init, psr) %>% 
  knitr::kable()
```

|  k\_stall|    k\_init|       psr|
|---------:|----------:|---------:|
|      0.02|  0.0039062|  0.003851|
|      0.10|  0.0039062|  0.003810|
|      0.50|  0.0039062|  0.003876|
|      0.02|  0.0078125|  0.007840|
|      0.10|  0.0078125|  0.007653|
|      0.50|  0.0078125|  0.007848|
|      0.02|  0.0156250|  0.015252|
|      0.10|  0.0156250|  0.015427|
|      0.50|  0.0156250|  0.015259|
|      0.02|  0.0312500|  0.019618|
|      0.10|  0.0312500|  0.030367|
|      0.50|  0.0312500|  0.030559|
|      0.02|  0.0625000|  0.019749|
|      0.10|  0.0625000|  0.059270|
|      0.50|  0.0625000|  0.059186|
|      0.02|  0.1250000|  0.019764|
|      0.10|  0.1250000|  0.093576|
|      0.50|  0.1250000|  0.111636|
|      0.02|  0.2500000|  0.019825|
|      0.10|  0.2500000|  0.093712|
|      0.50|  0.2500000|  0.199239|
|      0.02|  0.5000000|  0.019654|
|      0.10|  0.5000000|  0.093668|
|      0.50|  0.5000000|  0.326908|
|      0.02|  1.0000000|  0.019696|
|      0.10|  1.0000000|  0.093551|
|      0.50|  1.0000000|  0.367788|

Change in protein synthesis as a function of abort rate in SAT model
====================================================================

``` r
plot_data <- data %>% 
  mutate(k_stall = round(k_stall, 1)) %>%
  filter(k_stall == 0.1) %>% 
  filter(preterm_intact_rate < 0.1) %>% 
  mutate(preterm_intact_model = if_else(preterm_intact_rate == 0, "simple",
                                        preterm_intact_model)) %>% 
  filter(preterm_intact_model == "simple") %>% 
  mutate(model = forcats::fct_rev(model_names[preterm_intact_model]))

plot_data %>% 
  ggplot(aes(x = k_init, y = psr, color = as.factor(preterm_intact_rate), 
             shape = as.factor(preterm_intact_rate))) +
  geom_point(size = 1.5) + geom_line(size = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8.5, 2^0)) +
  scale_y_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))
                     ) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(19,17,15,18,16)) +
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       shape = "kabort (s-1)", color = "kabort (s-1)") +
  theme(legend.key.height = unit(0.2, "in")) +
                         
ggsave('figures/psr_sat_model_vary_abort_rate.pdf', width = 2.25, height = 1.45) 
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-9-1.png)

Source data for S1 Fig panel C
==============================

``` r
plot_data %>% 
  select(preterm_intact_rate, k_init, psr) %>% 
  knitr::kable()
```

|  preterm\_intact\_rate|    k\_init|       psr|
|----------------------:|----------:|---------:|
|                   0.00|  0.0039062|  0.003810|
|                   0.01|  0.0039062|  0.001859|
|                   0.02|  0.0039062|  0.000929|
|                   0.05|  0.0039062|  0.000077|
|                   0.00|  0.0078125|  0.007653|
|                   0.01|  0.0078125|  0.003700|
|                   0.02|  0.0078125|  0.001810|
|                   0.05|  0.0078125|  0.000168|
|                   0.00|  0.0156250|  0.015427|
|                   0.01|  0.0156250|  0.007207|
|                   0.02|  0.0156250|  0.003433|
|                   0.05|  0.0156250|  0.000387|
|                   0.00|  0.0312500|  0.030367|
|                   0.01|  0.0312500|  0.014169|
|                   0.02|  0.0312500|  0.006784|
|                   0.05|  0.0312500|  0.000708|
|                   0.00|  0.0625000|  0.059270|
|                   0.01|  0.0625000|  0.027856|
|                   0.02|  0.0625000|  0.013176|
|                   0.05|  0.0625000|  0.001432|
|                   0.00|  0.1250000|  0.093576|
|                   0.01|  0.1250000|  0.055900|
|                   0.02|  0.1250000|  0.025729|
|                   0.05|  0.1250000|  0.002797|
|                   0.00|  0.2500000|  0.093712|
|                   0.01|  0.2500000|  0.069528|
|                   0.02|  0.2500000|  0.046755|
|                   0.05|  0.2500000|  0.005243|
|                   0.00|  0.5000000|  0.093668|
|                   0.01|  0.5000000|  0.070689|
|                   0.02|  0.5000000|  0.051756|
|                   0.05|  0.5000000|  0.008220|
|                   0.00|  1.0000000|  0.093551|
|                   0.01|  1.0000000|  0.071091|
|                   0.02|  1.0000000|  0.052107|
|                   0.05|  1.0000000|  0.011020|

Change in protein synthesis as a function of stall elongation rate in SAT model
===============================================================================

``` r
plot_data <- data %>% 
  filter(preterm_intact_rate == 0.02) %>% 
  filter(preterm_intact_model == "simple") %>% 
  mutate(model = forcats::fct_rev(model_names[preterm_intact_model]))

plot_data %>% 
  ggplot(aes(x = k_init, y = psr, color = as.factor(k_stall), 
             shape = as.factor(k_stall))) +
  geom_point(size = 1.5) + geom_line(size = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8.5, 2^0)) +
  scale_y_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))
                     ) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(19,17,15,18,16)) +
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       shape = "kstall (s-1)", color = "kstall (s-1)") +
  theme(legend.key.height = unit(0.2, "in")) +
  geom_vline(aes(xintercept = k_stall, color = as.factor(k_stall)), 
             show.legend = F, linetype = "dotted")
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
ggsave('figures/psr_sat_model_vary_stall_rate.pdf', width = 2.25, height = 1.45) 
```

Source data for S1 Fig panel B
==============================

``` r
plot_data %>% 
  select(k_stall, k_init, psr) %>% 
  knitr::kable()
```

|  k\_stall|    k\_init|       psr|
|---------:|----------:|---------:|
|      0.02|  0.0039062|  0.000451|
|      0.10|  0.0039062|  0.000929|
|      0.50|  0.0039062|  0.001027|
|      0.02|  0.0078125|  0.000871|
|      0.10|  0.0078125|  0.001810|
|      0.50|  0.0078125|  0.002043|
|      0.02|  0.0156250|  0.001752|
|      0.10|  0.0156250|  0.003433|
|      0.50|  0.0156250|  0.004096|
|      0.02|  0.0312500|  0.003366|
|      0.10|  0.0312500|  0.006784|
|      0.50|  0.0312500|  0.007983|
|      0.02|  0.0625000|  0.006341|
|      0.10|  0.0625000|  0.013176|
|      0.50|  0.0625000|  0.015401|
|      0.02|  0.1250000|  0.007731|
|      0.10|  0.1250000|  0.025729|
|      0.50|  0.1250000|  0.031382|
|      0.02|  0.2500000|  0.007788|
|      0.10|  0.2500000|  0.046755|
|      0.50|  0.2500000|  0.061421|
|      0.02|  0.5000000|  0.007679|
|      0.10|  0.5000000|  0.051756|
|      0.50|  0.5000000|  0.102942|
|      0.02|  1.0000000|  0.007704|
|      0.10|  1.0000000|  0.052107|
|      0.50|  1.0000000|  0.166250|

Change in protein synthesis as a function of stall elongation rate in CSAT model
================================================================================

``` r
plot_data <- data %>% 
  filter(preterm_intact_rate == 1) %>% 
  filter(preterm_intact_model == "hit5") %>% 
  mutate(model = forcats::fct_rev(model_names[preterm_intact_model]))

plot_data %>% 
  ggplot(aes(x = k_init, y = psr, color = as.factor(k_stall), 
             shape = as.factor(k_stall))) +
  geom_point(size = 1.5) + geom_line(size = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     limits = c(2^-8.5, 2^0)) +
  scale_y_continuous(trans = "log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x, n = 4),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))
                     ) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(19,17,15,18,16)) +
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       shape = "kstall (s-1)", color = "kstall (s-1)") +
  theme(legend.key.height = unit(0.2, "in")) +
  geom_vline(aes(xintercept = k_stall, color = as.factor(k_stall)), 
             show.legend = F, linetype = "dotted")
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
ggsave('figures/psr_csat_model_vary_stall_rate.pdf', width = 2.25, height = 1.45) 
```

Source data for S1 Fig panel D
==============================

``` r
plot_data %>% 
  select(k_stall, k_init, psr) %>% 
  knitr::kable()
```

|  k\_stall|    k\_init|       psr|
|---------:|----------:|---------:|
|      0.02|  0.0039062|  0.003289|
|      0.10|  0.0039062|  0.003779|
|      0.50|  0.0039062|  0.003884|
|      0.02|  0.0078125|  0.005294|
|      0.10|  0.0078125|  0.007300|
|      0.50|  0.0078125|  0.007695|
|      0.02|  0.0156250|  0.007613|
|      0.10|  0.0156250|  0.013627|
|      0.50|  0.0156250|  0.014977|
|      0.02|  0.0312500|  0.007850|
|      0.10|  0.0312500|  0.023620|
|      0.50|  0.0312500|  0.028667|
|      0.02|  0.0625000|  0.004921|
|      0.10|  0.0625000|  0.035599|
|      0.50|  0.0625000|  0.052413|
|      0.02|  0.1250000|  0.001724|
|      0.10|  0.1250000|  0.041983|
|      0.50|  0.1250000|  0.088662|
|      0.02|  0.2500000|  0.000333|
|      0.10|  0.2500000|  0.035981|
|      0.50|  0.2500000|  0.131853|
|      0.02|  0.5000000|  0.000060|
|      0.10|  0.5000000|  0.024722|
|      0.50|  0.5000000|  0.166086|
|      0.02|  1.0000000|  0.000037|
|      0.10|  1.0000000|  0.019194|
|      0.50|  1.0000000|  0.181930|
