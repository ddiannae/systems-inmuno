library(GenomicDataCommons)
library(celldex)
library(SingleR)
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)

#Load ref
ref <- celldex::HumanPrimaryCellAtlasData()

# Load annotation file
annot_TCGA <- readRDS("annot_TCGA.RDS")

##### Download data ####

data_dir <- "/home/anakamura/.cache/GenomicDataCommons/"
dir.create(data_dir)
sample.sheet <- read.delim("gdc_sample_sheet.2023-06-27.tsv")
sample.sheet$Sample.ID <- gsub("-", "_", sample.sheet$Sample.ID)

fnames <- lapply(sample.sheet$File.ID, gdcdata)
# Would you like to create a GDC Cache directory at ~/.cache/GenomicDataCommons                                                                               
# 
# 1: Yes
# 2: No

#Select 1

#### Get Seurat object from sample i and normalization ####
for(i in 1:nrow(sample.sheet)){
  
  fileID_i <- sample.sheet[sample.sheet$File.ID == list.files(data_dir)[i],]$File.ID
  
  fileName_i <- sample.sheet[sample.sheet$File.ID == list.files(data_dir)[i],]$File.Name

  sample_i <- sample.sheet[sample.sheet$File.ID == list.files(data_dir)[i],]$Sample.ID
  
  file_dir <- paste0(data_dir, fileID_i, "/", fileName_i)
  
  destination_dir <- paste0("../../../datos/ot/anakamura/scRNA_npj_v3/", sample_i)
  dir.create(destination_dir)
  untar(file_dir, exdir = destination_dir)
  
  print(paste0("Get Seurat ", sample_i))
  
  FeatureBC_dir <- paste0(destination_dir, "/filtered_feature_bc_matrix")
  expression_matrix <- Read10X(data.dir = FeatureBC_dir, gene.column = 1)
  print(dim(expression_matrix))
  seurat_object = CreateSeuratObject(counts = expression_matrix, project = sample_i,
   min.cells = 100, min.features = 200)

##### Process Seurat i #####

print(paste0("Processing cells from ", sample_i))

#Create directories
dir.create(paste0(destination_dir, "/", sample_i, "_results"))
dir_res <- paste0(destination_dir, "/", sample_i, "_results")

#SCT normalization
seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", return.only.var.genes = FALSE)

#Get PCs
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# determine dimensionality of the data
elbow <- ElbowPlot(seurat_object) 
ggsave(filename =  paste0(dir_res, "/ElbowPlot_",sample_i, ".tiff"), plot = elbow, units="in",
       width=15, height=10, dpi =300, bg = "white")

seurat_object <- FindNeighbors(seurat_object, dims = 1:15)

# understanding resolution
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
#View(seurat_object@meta.data)

#Do UMAPs
seurat_object <- RunUMAP(seurat_object, dims = 1:15)
umap1 <- DimPlot(seurat_object, reduction = "umap")
ggsave(filename =  paste0(dir_res, "/ClustersUMAP_",sample_i, ".tiff"), plot = umap1, units="in", 
       width=20, height=10,  dpi =300)

#Idents(seurat_object) <- "SCT_snn_res.0.1"

#DEGs
seurat_object_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25,
                                     logfc.threshold = 0.1,  assay = "SCT")
seurat_object_markers <- seurat_object_markers %>% mutate(HGNC = annot_TCGA$gene_name[match(rownames(seurat_object_markers),
                                                                                      annot_TCGA$gene_id)])
write.table(seurat_object_markers, paste0(dir_res, "/markers_", sample_i,".tsv"), quote = FALSE, col.names = TRUE,
            row.names = TRUE)

#Get expression matrices
seurat_object.data <- GetAssayData(seurat_object, slot = "data", assay = "SCT")
seurat_object.cts <- GetAssayData(seurat_object, slot = "counts", assay = "SCT")
seurat_object.scale.data <- GetAssayData(seurat_object, slot = "scale.data", assay = "SCT")

original.rownames <- rownames(seurat_object.data)
rownames(seurat_object.data) <- annot_TCGA$gene_name[match(rownames(seurat_object.data), annot_TCGA$gene_id)]

#Annotate clusters by Main label
pred.main <- SingleR(seurat_object.data, 
                     ref = ref,
                     labels = ref$label.main, 
                     clusters = seurat_object$SCT_snn_res.0.5)
seurat_object@meta.data$singleR.labels.main <- pred.main$labels[match(seurat_object@meta.data$SCT_snn_res.0.5, rownames(pred.main))]
umap2 <- DimPlot(seurat_object, reduction = "umap", label = FALSE, 
                 group.by = "singleR.labels.main") 
ggsave(filename =  paste0(dir_res, "/CellsMainUMAP_",sample_i, ".tiff"), plot = umap2, units="in", 
       width=20, height=10,  dpi =300)

ht <- plotScoreHeatmap(pred.main)
ggsave(filename =  paste0(dir_res, "/MainScoreHeatmap_",sample_i, ".tiff"), plot = ht, units="in", 
       width=20, height=10,  dpi =300)

#Annotate clusters by Fine label
pred.fine <- SingleR(seurat_object.data, 
                     ref = ref,
                     labels = ref$label.fine, 
                     clusters = seurat_object$SCT_snn_res.0.5)
seurat_object@meta.data$singleR.labels.fine <- pred.fine$labels[match(seurat_object@meta.data$SCT_snn_res.0.5, rownames(pred.fine))]

umap3 <- DimPlot(seurat_object, reduction = "umap", label = FALSE, 
                 group.by ="singleR.labels.fine")
ggsave(filename =  paste0(dir_res, "/CellsFineUMAP_",sample_i, ".tiff"), plot = umap3, units="in", 
       width=20, height=10,  dpi =300)

ht2 <- plotScoreHeatmap(pred.fine)
ggsave(filename =  paste0(dir_res, "/FineScoreHeatmap_",sample_i, ".tiff"), plot = ht2, units="in", 
       width=20, height=10,  dpi =300)

#Get matrix per cluster

dir.create(paste0(dir_res, "/expression_matrices"))

for(j in 1:length(unique(seurat_object@meta.data$SCT_snn_res.0.5))){
  
  print(paste0("Saving matrix from ", sample_i, " cluster ", 
        as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j])))
  
  assign("matrix.cts",
         seurat_object.cts[,which(colnames(seurat_object.cts) %in% 
                                 rownames(seurat_object@meta.data[seurat_object@meta.data$SCT_snn_res.0.5 == as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),
                                 ]))
         ])
  
  matrix.cts <- as.data.frame(matrix.cts)
  
  assign("matrix.data",
         seurat_object.data[,which(colnames(seurat_object.data) %in% 
                                    rownames(seurat_object@meta.data[seurat_object@meta.data$SCT_snn_res.0.5 == as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),
                                    ]))
         ])
  
  matrix.data <- as.data.frame(matrix.data)
  
  assign("matrix.scale.data",
         seurat_object.scale.data[,which(colnames(seurat_object.scale.data) %in% 
                                    rownames(seurat_object@meta.data[seurat_object@meta.data$SCT_snn_res.0.5 == as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),
                                    ]))
         ])
  
  matrix.scale.data <- as.data.frame(matrix.scale.data)
  
  #Prep rownames to Ensembl
  original.rownames.NoVersion <- gsub("\\.[0-9]*", "", original.rownames)

  rownames(matrix.cts) <- original.rownames.NoVersion
  rownames(matrix.data) <- original.rownames.NoVersion
  rownames(matrix.scale.data) <- original.rownames.NoVersion
  
  #Agregate Column Gene as first column
  matrix.cts <- matrix.cts %>% mutate(Gene = rownames(matrix.cts)) %>% dplyr::relocate(Gene)
  matrix.data <- matrix.data %>% mutate(Gene = rownames(matrix.data)) %>% dplyr::relocate(Gene)
  matrix.scale.data <- matrix.scale.data %>% mutate(Gene = rownames(matrix.scale.data)) %>% dplyr::relocate(Gene)
  
  write.table(matrix.cts,
              file = paste0(dir_res, "/expression_matrices/Cls_",as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),"_Counts_",sample_i,".tsv"),
              quote = FALSE,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE)
  
  write.table(matrix.data,
              file = paste0(dir_res, "/expression_matrices/Cls_",as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),"_Data_",sample_i,".tsv"),
              quote = FALSE,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE)
  
  write.table(matrix.scale.data,
              file = paste0(dir_res, "/expression_matrices/Cls_",as.numeric(levels(seurat_object@meta.data$SCT_snn_res.0.5)[j]),"_ScaleData_",sample_i,".tsv"),
              quote = FALSE,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE)

  
}

#Save metadata
write.csv(seurat_object@meta.data, file = paste0(dir_res, "/metadata_", sample_i, ".csv"), 
          quote = FALSE)
#Save Seurat obj with processed data
saveRDS(seurat_object, paste0(dir_res, "/Processed_Seurat_",sample_i,".RDS"))
}
