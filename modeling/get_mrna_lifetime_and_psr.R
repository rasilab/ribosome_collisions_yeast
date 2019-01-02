#!/fh/fast/subramaniam_a/user/rasi/lib/R-devel/bin/Rscript

#'Read in reaction firings and calculate statistics of various reaction firing
#'frequencies.
#'Currently looking at mRNA lifetime, protein synthesis rate, collision rate,
#'premature termination rate
#'@author Arvind R. Subramniam
#'@date 29 Aug 2018

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly=TRUE)
# which simulation to look at
current_sim <- args[1]
# this is the folder from which the rxns.tsv file is read and the output written to
output_folder <- args[2]

filter_fread <- function (file) {
  sim_id = as.integer(str_extract(file, "(?<=tasep_)[:digit:]+"))
  suppressMessages(read_tsv(file)) %>% 
    filter(str_detect(rxn, "transcription|exonucleolysis_5end_0|^term|^collision|^preterm"))
}

all_data <- list.files(output_folder, pattern = "rxns.tsv", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(sim_id = as.integer(str_extract(file, "(?<=tasep_)[:digit:]+"))) %>%
  filter(sim_id == as.integer(current_sim)) %>% 
  # take out only the reactions of interest to keep dataframe small
  mutate(data = map(file, filter_fread)) %>%
  select(-file, -sno) %>%
  unnest()

lifetime_data <- all_data %>%
  # lifetime is defined as the time from transcription to iniitation of exonu-
  # cleolytic 5' degradation
  filter(str_detect(rxn, "transcription|exonucleolysis_5end_0")) %>%
  group_by(sim_id, mol_id) %>%
  summarize(lifetime = max(time) - min(time), num_events = n()) %>%
  ungroup() %>%
  # this will be 0 if only transcription or exonucleolytic decay was observed
  # also there should be only one transcription and exonucleolysis event
  filter(lifetime > 0 & num_events == 2) %>%
  # get the mean and sd of mRNA lifetime in each simulation along with how many
  # mRNAs were observed
  group_by(sim_id) %>%
  summarize(mean_lifetime = as.integer(mean(lifetime)),
            sd_lifetime = as.integer(sd(lifetime)), n_mrna = n()) %>%
  ungroup() %>%
  write_tsv(paste0(output_folder, "/tasep_", current_sim, "_mrna_lifetime_stats.tsv.gz"))

psr_data <- all_data %>% 
  group_by(sim_id, mol_id) %>% 
  # first calculate the maximum time in the simulation for each mRNA molecule
  mutate(time = max(time)) %>% 
  ungroup() %>% 
  # extract all normal termination events
  filter(str_detect(rxn, "^term")) %>% 
  # calculate proteins per mRNA molecule
  group_by(sim_id, mol_id) %>% 
  summarize(n_protein = n(), time = max(time)) %>% 
  ungroup() %>% 
  # calculate mean and sd of proteins per mRNA and average protein synthesis
  # rate across the whole simulation
  group_by(sim_id) %>% 
  summarize(mean_p_per_m = as.integer(mean(n_protein)), 
            sd_p_per_m = as.integer(sd(n_protein)), 
            total_p = sum(n_protein),
            total_time = as.integer(max(time)),
            psr = round(sum(n_protein)/max(time), 6)) %>% 
  ungroup() %>% 
  write_tsv(paste0(output_folder, "/tasep_", current_sim, "_psr_stats.tsv.gz"))

all_data %>% 
  group_by(sim_id, mol_id) %>% 
  # first calculate the maximum time in the simulation for each mRNA molecule
  mutate(time = max(time)) %>% 
  ungroup() %>% 
  # extract all normal termination events
  filter(str_detect(rxn, "^collision")) %>% 
  # calculate collisions per mRNA molecule
  group_by(sim_id, mol_id) %>% 
  summarize(n_collision = n(), time = max(time)) %>% 
  ungroup() %>% 
  # calculate mean and sd of collisions per mRNA and average collision synthesis
  # rate across the whole simulation
  group_by(sim_id) %>% 
  summarize(mean_p_per_m = as.integer(mean(n_collision)), 
            sd_p_per_m = as.integer(sd(n_collision)), 
            total_collision = sum(n_collision),
            total_time = as.integer(max(time)),
            collision_freq = round(sum(n_collision)/max(time), 6)) %>% 
  ungroup() %>% 
  write_tsv(paste0(output_folder, "/tasep_", current_sim, "_collision_stats.tsv.gz"))

all_data %>% 
  group_by(sim_id, mol_id) %>% 
  # first calculate the maximum time in the simulation for each mRNA molecule
  mutate(time = max(time)) %>% 
  ungroup() %>% 
  # extract all premature termination events
  filter(str_detect(rxn, "^preterm")) %>% 
  # extract premature_term_type
  mutate(rxn_type = str_extract(rxn, ".+(?=_\\d+)")) %>% 
  # calculate num of premature termination events of each type per mRNA molecule
  group_by(sim_id, mol_id, rxn_type) %>% 
  summarize(n_abort = n(), time = max(time)) %>% 
  ungroup() %>% 
  # calculate mean and sd of premature terms per mRNA and average premature term
  # rate across the whole simulation
  group_by(sim_id, rxn_type) %>% 
  summarize(mean_abort_per_m = as.integer(mean(n_abort)), 
            sd_abort_per_m = as.integer(sd(n_abort)), 
            total_abort = sum(n_abort),
            total_time = as.integer(max(time)),
            abort_freq = round(sum(n_abort)/max(time), 6)) %>% 
  ungroup() %>% 
  write_tsv(paste0(output_folder, "/tasep_", current_sim, "_abort_stats.tsv.gz"))
