library(readr)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

cat("Reading files\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Three arguments must be provided", call.=FALSE)
  
}

MEM <- args[1]
EXPR <- args[2]
OUT <- args[3]

membership <- read_tsv(MEM)
gene_universe <- read_tsv(EXPR, col_types = cols_only(gene_id = col_character()) )

gene_universe <- gene_universe %>% 
  separate(col = gene_id, into = c("ensembl", "version"), sep = "\\.") %>%
  pull(ensembl)

membership <- membership %>% 
  separate(col = ensembl, into = c("ensembl", "version"), sep = "\\.") %>%
  select(-version)

all_enrichments <- lapply(X = unique(membership$community),
                          FUN = function(com){
                            
                            cat("Working with community: ", com, "\n")
                            
                            gene_list <- membership %>% filter(community == com) %>%
                              pull(ensembl)
                            
                            if(length(gene_list) >= 5) {
                              
                              terms <- enrichGO(gene          = gene_list,
                                                universe      = gene_universe,
                                                OrgDb         = org.Hs.eg.db,
                                                keyType       = "ENSEMBL",
                                                ont           = "BP",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.005,
                                                qvalueCutoff  = 0.01,
                                                minGSSize     = 10,
                                                readable      = FALSE)
                              
                              tryCatch({
                                bp <- pairwise_termsim(terms)
                                bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
                                simple_results <- bp2@result
                                if(nrow(simple_results) > 0) {
                                  simple_results$commun <- com 
                                }
                                return(simple_results)
                              }, error = function(cond) {
                                return(NULL)
                              })
                              
                            }
                            return(NULL)
                          })

bind_rows(all_enrichments) %>% 
  janitor::clean_names() %>%
  write_tsv(OUT)


