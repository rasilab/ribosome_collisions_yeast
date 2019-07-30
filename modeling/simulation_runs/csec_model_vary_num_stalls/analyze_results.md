Analyze simulation results
================
rasi
30 July, 2019

-   [Load libraries](#load-libraries)
-   [Read protein count data](#read-protein-count-data)
-   [Read simulation parameters](#read-simulation-parameters)
-   [Combine all data into a single table](#combine-all-data-into-a-single-table)
-   [PSR as a function of number of stalls for supplementary figure](#psr-as-a-function-of-number-of-stalls-for-supplementary-figure)
-   [Source data for S3 Fig panel B](#source-data-for-s3-fig-panel-b)

Load libraries
--------------

``` r
library(tidyverse)
library(rasilabRtemplates)
# disable scientific notation
options(scipen=999)

cleave_model_names <- c(
  "hit5" = "CSEC",
  "simple" = "SEC",
  "trafficjam" = "TJ"
)
```

Read protein count data
=======================

``` r
psr_data <- read_tsv("tables/psr_stats.tsv") %>% 
  print()
```

    ## # A tibble: 90 x 6
    ##    sim_id mean_p_per_m sd_p_per_m total_p total_time     psr
    ##     <int>        <int>      <int>   <int>      <int>   <dbl>
    ##  1      0            3          2    3315     999930 0.00332
    ##  2      1            6          4    6650     999561 0.00665
    ##  3     10            6          4    6876     999732 0.00688
    ##  4     11           13          7   13329     999929 0.0133 
    ##  5     12           26         13   27553     999995 0.0276 
    ##  6     13           49         25   47056     997522 0.0472 
    ##  7     14           40         29   39676     999997 0.0397 
    ##  8     15           19         11   19931     998784 0.0200 
    ##  9     16           18         12   17731     999461 0.0177 
    ## 10     17           20         12   20974     998908 0.0210 
    ## # ... with 80 more rows

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
  # retain only parameters that are varied, the others are for checking
  group_by(parameter) %>% 
  mutate(vary = if_else(length(unique(value)) > 1, T, F)) %>% 
  ungroup() %>% 
  filter(vary == T) %>% 
  select(-vary) %>% 
  spread(parameter, value) %>% 
  print()
```

    ## # A tibble: 90 x 12
    ##    sim_id k_elong_stall_1 k_elong_stall_2 k_elong_stall_3 k_elong_stall_4
    ##     <int>           <dbl>           <dbl>           <dbl>           <dbl>
    ##  1      0             0.1            NA                NA              NA
    ##  2      1             0.1            NA                NA              NA
    ##  3      2             0.1            NA                NA              NA
    ##  4      3             0.1            NA                NA              NA
    ##  5      4             0.1            NA                NA              NA
    ##  6      5             0.1            NA                NA              NA
    ##  7      6             0.1            NA                NA              NA
    ##  8      7             0.1            NA                NA              NA
    ##  9      8             0.1            NA                NA              NA
    ## 10      9             0.2             0.2              NA              NA
    ## # ... with 80 more rows, and 7 more variables: k_elong_stall_5 <dbl>,
    ## #   k_elong_stall_6 <dbl>, k_elong_stall_7 <dbl>, k_elong_stall_8 <dbl>,
    ## #   k_elong_stall_9 <dbl>, k_init <dbl>, n_stall <dbl>

Combine all data into a single table
====================================

``` r
data <- annotations %>% 
  left_join(psr_data, by = "sim_id") %>% 
  filter(n_stall <= 6) %>%
  mutate(n_stall = as.factor(n_stall)) %>% 
  print()
```

    ## # A tibble: 54 x 17
    ##    sim_id k_elong_stall_1 k_elong_stall_2 k_elong_stall_3 k_elong_stall_4
    ##     <int>           <dbl>           <dbl>           <dbl>           <dbl>
    ##  1      0             0.1            NA                NA              NA
    ##  2      1             0.1            NA                NA              NA
    ##  3      2             0.1            NA                NA              NA
    ##  4      3             0.1            NA                NA              NA
    ##  5      4             0.1            NA                NA              NA
    ##  6      5             0.1            NA                NA              NA
    ##  7      6             0.1            NA                NA              NA
    ##  8      7             0.1            NA                NA              NA
    ##  9      8             0.1            NA                NA              NA
    ## 10      9             0.2             0.2              NA              NA
    ## # ... with 44 more rows, and 12 more variables: k_elong_stall_5 <dbl>,
    ## #   k_elong_stall_6 <dbl>, k_elong_stall_7 <dbl>, k_elong_stall_8 <dbl>,
    ## #   k_elong_stall_9 <dbl>, k_init <dbl>, n_stall <fct>,
    ## #   mean_p_per_m <int>, sd_p_per_m <int>, total_p <int>, total_time <int>,
    ## #   psr <dbl>

PSR as a function of number of stalls for supplementary figure
==============================================================

``` r
plot_data <- data

plot_data %>%
  ggplot(aes(x = k_init, y = psr, color = n_stall, group = n_stall)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = "log2",
                     labels = scales::trans_format("log2", scales::math_format(2^.x)),
                     breaks = 2^(seq(-8,0,2))) +
  scale_y_continuous(limits = c(0, 0.06)) +
  viridis::scale_color_viridis(discrete = T, end = 0.9) +
  labs(x = "initiation rate (s-1)", y = "protein synthesis rate (s-1)",
       color = "number\nof stalls", shape = "") +
  theme(legend.key.height = unit(0.2, "in")) +
  guides(color = guide_legend(
                 keywidth=0.1,
                 keyheight=0.15,
                 default.unit="inch")
      )
```

![](analyze_results_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
ggsave("figures/psr_vs_initiation_rate_vary_n_stalls.pdf")
```

Source data for S3 Fig panel B
==============================

``` r
plot_data %>% 
  select(n_stall, k_init, psr) %>% 
  knitr::kable()
```

| n\_stall |    k\_init|       psr|
|:---------|----------:|---------:|
| 1        |  0.0039062|  0.003315|
| 1        |  0.0078125|  0.006653|
| 1        |  0.0156250|  0.013309|
| 1        |  0.0312500|  0.026320|
| 1        |  0.0625000|  0.044162|
| 1        |  0.1250000|  0.038605|
| 1        |  0.2500000|  0.018621|
| 1        |  0.5000000|  0.020688|
| 1        |  1.0000000|  0.020902|
| 2        |  0.0039062|  0.003454|
| 2        |  0.0078125|  0.006878|
| 2        |  0.0156250|  0.013330|
| 2        |  0.0312500|  0.027553|
| 2        |  0.0625000|  0.047173|
| 2        |  0.1250000|  0.039676|
| 2        |  0.2500000|  0.019955|
| 2        |  0.5000000|  0.017741|
| 2        |  1.0000000|  0.020997|
| 3        |  0.0039062|  0.003476|
| 3        |  0.0078125|  0.006950|
| 3        |  0.0156250|  0.013578|
| 3        |  0.0312500|  0.026602|
| 3        |  0.0625000|  0.049207|
| 3        |  0.1250000|  0.042888|
| 3        |  0.2500000|  0.018867|
| 3        |  0.5000000|  0.018156|
| 3        |  1.0000000|  0.019478|
| 4        |  0.0039062|  0.003397|
| 4        |  0.0078125|  0.006864|
| 4        |  0.0156250|  0.013187|
| 4        |  0.0312500|  0.026930|
| 4        |  0.0625000|  0.048940|
| 4        |  0.1250000|  0.044049|
| 4        |  0.2500000|  0.021119|
| 4        |  0.5000000|  0.018760|
| 4        |  1.0000000|  0.019611|
| 5        |  0.0039062|  0.003510|
| 5        |  0.0078125|  0.006771|
| 5        |  0.0156250|  0.013581|
| 5        |  0.0312500|  0.026761|
| 5        |  0.0625000|  0.049855|
| 5        |  0.1250000|  0.047190|
| 5        |  0.2500000|  0.020213|
| 5        |  0.5000000|  0.019033|
| 5        |  1.0000000|  0.021675|
| 6        |  0.0039062|  0.003503|
| 6        |  0.0078125|  0.006808|
| 6        |  0.0156250|  0.013802|
| 6        |  0.0312500|  0.027153|
| 6        |  0.0625000|  0.050775|
| 6        |  0.1250000|  0.045152|
| 6        |  0.2500000|  0.020907|
| 6        |  0.5000000|  0.019538|
| 6        |  1.0000000|  0.020357|
