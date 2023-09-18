library(readr)
library(dplyr)
library(janitor)

args <- commandArgs(trailingOnly=TRUE)
COUNTSDIR <- args[1]
ARACNEDIR <-  args[2]

getAracneMatrix <- function(counts, samples) {
  
  counts %>% 
    select(gene_id, all_of(samples))
}

normalized_counts <- read_tsv(paste0(COUNTSDIR, "/", "normalized_counts.tsv"))
vst_counts <- read_tsv(paste0(COUNTSDIR, "/", "vst_counts.tsv"))
cases_data <- read_tsv(paste0(COUNTSDIR, "/", "cases_data.tsv"))
controls_data <- read_tsv(paste0(COUNTSDIR, "/", "controls_data.tsv"))

getAracneMatrix(normalized_counts, cases_data$id) %>%
  write_tsv(paste0(ARACNEDIR, "/normalized_cases.tsv"))

getAracneMatrix(normalized_counts, controls_data$id) %>%
  write_tsv(paste0(ARACNEDIR, "/normalized_controls.tsv"))

getAracneMatrix(vst_counts, cases_data$id) %>%
  write_tsv(paste0(ARACNEDIR, "/vst_cases.tsv"))

getAracneMatrix(vst_counts, controls_data$id) %>%
  write_tsv(paste0(ARACNEDIR, "/vst_controls.tsv"))
