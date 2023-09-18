library(vroom)
library(dplyr)
library(tidyr)
library(janitor)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("Four arguments must be provided", call.=FALSE)
}

ADJ <- args[1]
CUTOFF <-  as.numeric(args[2])
ANNOT <-  args[3]
OUT <-  args[4]
annot <- vroom(ANNOT) %>%
  clean_names()

cat("Reading matrix\n")
MIvals <- vroom(ADJ)
genes <- colnames(MIvals)

cat("MI matrix with ", nrow(MIvals), " rows and ", ncol(MIvals), " columns loaded \n")
MIvals$source <- genes

cat("Getting data frame\n")
MIvals <- MIvals %>% 
  pivot_longer(cols = starts_with("ENSG"), names_to = "target",
               values_to = "mi", values_drop_na = TRUE) %>%
  filter(target != source) %>%
  arrange(desc(mi)) %>% 
  mutate(row_num = row_number()) %>%
  filter(row_num <= CUTOFF)

cat("Annotating interactions\n")
annot <- annot %>% 
  select(gene_id = gene_stable_id_version, chr = chromosome_scaffold_name, 
         gene_type = gene_type, symbol = hgnc_symbol)
colnames(annot) <-  c("source", "source_chr", "source_type", "source_symbol")
MIvals <- merge(annot, MIvals)

colnames(annot) <-  c("target", "target_chr", "target_type", "target_symbol")
MIvals <- merge(annot, MIvals)
MIvals <- MIvals %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                            interaction_type = if_else(inter == T, "Inter", "Intra"))

targets <- MIvals %>% select(target, target_chr, target_type, target_symbol)
sources <- MIvals %>% select(source, source_chr, source_type, source_symbol)
colnames(targets) <- c("ensembl", "chr", "type",  "symbol")
colnames(sources)  <-  c("ensembl", "chr", "type",  "symbol")
vertices <- bind_rows(targets, sources)
vertices <- vertices[!duplicated(vertices$ensembl), ]

MIvals <- MIvals %>% 
  select(source, target, mi, row_num, interaction_type) %>%
  arrange(row_num) 

cat("Saving files\n")
vroom_write(vertices, file = paste0(OUT, "_vertices.tsv"))
vroom_write(MIvals, file = paste0(OUT, "_interactions.tsv"))