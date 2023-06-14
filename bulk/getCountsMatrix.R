library(readr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Two arguments must be provided", call.=FALSE)
}

DIR <- args[1]
TYPE <- args[2]
cat("Reading files in", paste0(DIR, "/", TYPE), "\n")

files_to_read <- list.files(path = paste0(DIR, "/", TYPE), 
                            pattern = "\\.rna_seq.augmented_star_gene_counts.tsv$", full.names = T, recursive = T)

cat("Found", length(files_to_read), "files\n")

targets <- tibble(id = paste(TYPE, 1:length(files_to_read), sep = "_"),
                      file = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 1))),
                      file_id = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 2)[1])),
                      group = TYPE)

cat("Reading files\n")

all_data <- lapply(files_to_read, function(file) {
  id <- targets %>%
    filter(file == tail(strsplit(!!file, "/")[[1]], n = 1)) %>%
    pull(id)
  data <- read_tsv(file, comment = "#") %>% 
    select(gene_id, gene_name, gene_type, !!id := unstranded) %>% 
    filter(!startsWith(gene_id, "N_"))
  
  return(data)
})                            

cat("Writing data\n")
all_data <- Reduce(inner_join, all_data) 
all_data %>%
  write_tsv(paste0(DIR, "/", TYPE, "_counts.tsv"))

targets %>% 
  write_tsv(paste0(DIR, "/", TYPE, "_data.tsv"))
  


