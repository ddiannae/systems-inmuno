library(celldex)
library(SingleR)
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)

#Set working directory in the same as the one used to download the data
setwd("~/scRNAseq_data/")

#Load ref
ref <- celldex::HumanPrimaryCellAtlasData()

# Load annotation file
annot_TCGA <- readRDS("annot_TCGA.RDS")

#Load Seurat object
files <- list.files(pattern = "seurat.*RDS")

#Get sample name
n <- gsub("^seurat_", "", files)
n <- gsub("\\.RDS$", "", n)
n

# Iterate over every element on the list
for(i in 1:length(n)) {
  
  print(paste0("Processing cells from ", n[i]))
  
  #Create directories
  dir.create(paste0("results/",n[i], "_results"))
  dir_res <- paste0("results/",n[i], "_results")
  
  #Load Seurat
  seurat.sc1 <- readRDS(files[i])
  
  seurat.sc1 <- SCTransform(seurat.sc1)
  
  seurat.sc1 <- RunPCA(seurat.sc1, features = VariableFeatures(object = seurat.sc1))
  
  # determine dimensionality of the data
  elbow <- ElbowPlot(seurat.sc1) 
  ggsave(filename =  paste0(dir_res, "/ElbowPlot_",n[i], ".tiff"), plot = elbow, units="in",
         width=15, height=10, dpi =300, bg = "white")
  
  seurat.sc1 <- FindNeighbors(seurat.sc1, dims = 1:15)
  
  # understanding resolution
  seurat.sc1 <- FindClusters(seurat.sc1, resolution = 0.5)
  #View(seurat.sc1@meta.data)
  
  seurat.sc1 <- RunUMAP(seurat.sc1, dims = 1:15)
  umap1 <- DimPlot(seurat.sc1, reduction = "umap")
  ggsave(filename =  paste0(dir_res, "/ClustersUMAP_",n[i], ".tiff"), plot = umap1, units="in", 
         width=20, height=10,  dpi =300)
  
  #Idents(seurat.sc1) <- "SCT_snn_res.0.1"
  
  seurat.sc1_markers <- FindAllMarkers(seurat.sc1, only.pos = TRUE, min.pct = 0.25,
                                       logfc.threshold = 0.1)
  seurat.sc1_markers <- seurat.sc1_markers %>% mutate(HGNC = annot_TCGA$gene_name[match(rownames(seurat.sc1_markers),
                                                                                        annot_TCGA$gene_id)])
  write.table(seurat.sc1_markers, paste0(dir_res, "/markers_", n[i],".tsv"), quote = FALSE, col.names = TRUE,
              row.names = TRUE)
  
  #Annotate clusters
  seurat.sc1.cts <- GetAssayData(seurat.sc1)
  original.rownames <- rownames(seurat.sc1.cts)
  rownames(seurat.sc1.cts) <- annot_TCGA$gene_name[match(rownames(seurat.sc1.cts), annot_TCGA$gene_id)]
  
  #Main label
  pred.main <- SingleR(seurat.sc1.cts, 
                       ref = ref,
                       labels = ref$label.main, 
                       clusters = seurat.sc1$SCT_snn_res.0.5)
  seurat.sc1@meta.data$singleR.labels.main <- pred.main$labels[match(seurat.sc1@meta.data$SCT_snn_res.0.5, rownames(pred.main))]
  umap2 <- DimPlot(seurat.sc1, reduction = "umap", label = FALSE, 
                   group.by = "singleR.labels.main") 
  ggsave(filename =  paste0(dir_res, "/CellsMainUMAP_",n[i], ".tiff"), plot = umap2, units="in", 
         width=20, height=10,  dpi =300)
  
  ht <- plotScoreHeatmap(pred.main)
  ggsave(filename =  paste0(dir_res, "/MainScoreHeatmap_",n[i], ".tiff"), plot = ht, units="in", 
         width=20, height=10,  dpi =300)
  
  #Fine label
  pred.fine <- SingleR(seurat.sc1.cts, 
                       ref = ref,
                       labels = ref$label.fine, 
                       clusters = seurat.sc1$SCT_snn_res.0.5)
  seurat.sc1@meta.data$singleR.labels.fine <- pred.fine$labels[match(seurat.sc1@meta.data$SCT_snn_res.0.5, rownames(pred.fine))]
  
  umap3 <- DimPlot(seurat.sc1, reduction = "umap", label = FALSE, 
                   group.by ="singleR.labels.fine")
  ggsave(filename =  paste0(dir_res, "/CellsFineUMAP_",n[i], ".tiff"), plot = umap3, units="in", 
         width=20, height=10,  dpi =300)
  
  ht2 <- plotScoreHeatmap(pred.fine)
  ggsave(filename =  paste0(dir_res, "/FineScoreHeatmap_",n[i], ".tiff"), plot = ht2, units="in", 
         width=20, height=10,  dpi =300)
  

  #Get matrix per cluster
  
  dir.create(paste0(dir_res, "/expression_matrices"))
  
  for(j in 1:length(unique(seurat.sc1@meta.data$SCT_snn_res.0.5))){
    
    assign("matrix",
           seurat.sc1.cts[,which(colnames(seurat.sc1.cts) %in% 
          rownames(seurat.sc1@meta.data[seurat.sc1@meta.data$SCT_snn_res.0.5 == as.numeric(levels(seurat.sc1@meta.data$SCT_snn_res.0.5)[j]),
                                   ]))
           ])
    
    matrix <- as.data.frame(matrix)
    
    #Prep rownames to Ensembl
    original.rownames.NoVersion <- gsub("\\.[0-9]*", "", original.rownames)
    rownames(matrix) <- original.rownames.NoVersion
    
    #Filter
    threshold <- round(dim(matrix)[2]/2)
    ridx <- rowSums(matrix == 0) <= threshold
    matrix1 <- matrix[ridx, ]
    
    #Agregate Column Gene as first column
    matrix <- matrix %>% mutate(Gene = rownames(matrix)) %>% dplyr::relocate(Gene)
    matrix1 <- matrix1 %>% mutate(Gene = rownames(matrix1)) %>% dplyr::relocate(Gene)
    #Change Gene column to Ensembl
    # matrix$Gene <- annot_TCGA$Ensembl_noVers[match(rownames(matrix), annot_TCGA$gene_name)]
    # matrix1$Gene <- annot_TCGA$Ensembl_noVers[match(rownames(matrix1), annot_TCGA$gene_name)]
    
    write.table(matrix,
                file = paste0(dir_res, "/expression_matrices/Cls_",unique(seurat.sc1@meta.data$SCT_snn_res.0.5)[j],"_cts_",n[i],".tsv"),
                quote = FALSE,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE)
    
    write.table(matrix1,
                file = paste0(dir_res, "/expression_matrices/Cls_",unique(seurat.sc1@meta.data$SCT_snn_res.0.5)[j],"_FilteredCts_",n[i],".tsv"),
                quote = FALSE,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE)
    
  }
  
  
  #Save metadata
  write.csv(seurat.sc1@meta.data, file = paste0(dir_res, "/metadata_", n[i], ".csv"), 
            quote = FALSE)
  #Save Seurat obj with processed data
  saveRDS(seurat.sc1, paste0(dir_res, "/Processed_Seurat_",n[i],".RDS"))
}
