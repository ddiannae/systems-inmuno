library(readr)
library(dplyr)
library(janitor)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Two arguments must be provided", call.=FALSE)
}

DIR <- args[1]
TYPE <- args[2]
ANNOT <- args[3]

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

all_data <- Reduce(inner_join, all_data) 

all_data <- all_data %>%
  filter(gene_type %in% c("protein_coding")) %>%
  distinct(gene_id, .keep_all = T) %>%
  distinct(gene_name, .keep_all = T) 

annot <- read_tsv(ANNOT) %>%
  clean_names() %>%
  filter(chromosome_scaffold_name %in% c(as.character(1:22), "X", "Y")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  select(gene_id = gene_stable_id_version, gene_name = hgnc_symbol, everything())

all_data <- all_data %>%
  semi_join(annot, by = c("gene_id", "gene_name"))

head(cat("Writing data\n"))

all_data %>%
  write_tsv(paste0(DIR, "/", TYPE, "_counts.tsv"))

targets %>% 
  write_tsv(paste0(DIR, "/", TYPE, "_data.tsv"))
  