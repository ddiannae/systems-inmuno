library(readr)
library(tibble)
library(dplyr)
library(DESeq2)
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)
COUNTSDIR <- args[1]
PLOTSDIR <- args[2]

# Reading files
cases_counts <- read_tsv(paste0(COUNTSDIR, "/", "cases_counts.tsv"))
controls_counts <- read_tsv(paste0(COUNTSDIR, "/", "controls_counts.tsv"))
cases_data <- read_tsv(paste0(COUNTSDIR, "/", "cases_data.tsv"))
controls_data <- read_tsv(paste0(COUNTSDIR, "/", "controls_data.tsv"))

# Binding data tibbles
all_data <- bind_rows(controls_data, cases_data)

# Joining counts
all_counts <- controls_counts %>% 
  inner_join(cases_counts) %>% 
  select(gene_id, !!all_data$id) %>%
  column_to_rownames(var = "gene_id") %>% 
  as.matrix()

# Filtering low expression genes
exp_genes <- apply(all_counts, 1, function(x) mean(x) > 10)
egtable <- table(exp_genes)
cat("There are", egtable[[1]], "genes with mean expression count < 10 and ", egtable[[2]], "with mean count > 10 \n")

non_zero_genes <- apply(all_counts, 1, function(x) sum(x == 0) < ncol(all_counts)/2)
egtable <- table(non_zero_genes)
cat("There are", egtable[[1]], "genes with 0 values in more than 50% samples", egtable[[2]], 
    " genes with 0 values in less than 50% samples \n")

all_counts <- all_counts[non_zero_genes & exp_genes,]

# Normalization with DESeq2
dds <- DESeqDataSetFromMatrix(countData = all_counts,
                              colData = column_to_rownames(all_data, var = "id"), 
                              design = ~group)
dds <- estimateSizeFactors(dds)

# Saving normalized counts
counts(dds, normalized=TRUE) %>% 
  as_tibble( rownames = "gene_id") %>%
  write_tsv(paste0(COUNTSDIR, "/normalized_counts.tsv"))

# Saving vst counts
vst <- vst(dds, blind=TRUE)
vst_counts <- assay(vst)
vst_counts %>% 
  as_tibble( rownames = "gene_id") %>%
  write_tsv(paste0(COUNTSDIR, "/vst_counts.tsv"))

# Saving PCA
png(paste(PLOTSDIR, "vst_pca.png", sep="/"), width = 800)
plotPCA(vst, intgroup="group")
dev.off()

vst_cor <- cor(vst_counts)  

# Saving heatmap
png(paste(PLOTSDIR, "vst_heatmap.png", sep="/"), width = 1200, height = 1200)
pheatmap(vst_cor)
dev.off()

# Saving DeSeq2 object
saveRDS(dds, paste0(COUNTSDIR, "/dds.rds"))

        