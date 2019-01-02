#!/fh/fast/subramaniam_a/user/rasi/lib/R-devel/bin/Rscript

#' Combine mRNA lifetime and protein data for all simulations
#' @author Arvind R. Subramaniam 
#' @date Thu Aug 30 07:47:19 2018 


library(tidyverse)

list.files("output/", pattern = "mrna_lifetime_stats.tsv.gz", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, read_tsv)) %>% 
  select(-sno, -file) %>% 
  unnest() %>% 
  write_tsv("tables/mrna_lifetime_stats.tsv")

list.files("output/", pattern = "_psr_stats.tsv.gz", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, read_tsv)) %>% 
  select(-sno, -file) %>% 
  unnest() %>% 
  write_tsv("tables/psr_stats.tsv")

list.files("output/", pattern = "_collision_stats.tsv.gz", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, read_tsv)) %>% 
  select(-sno, -file) %>% 
  unnest() %>% 
  write_tsv("tables/collision_stats.tsv")

list.files("output/", pattern = "_abort_stats.tsv.gz", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, read_tsv)) %>% 
  # get rid of empty dataframes
  mutate(nrows = map_int(data, nrow)) %>% 
  filter(nrows > 0) %>% 
  select(-sno, -file, nrows) %>%
  unnest() %>%
  write_tsv("tables/abort_stats.tsv")
