#Load Sample Sheet
sample.sheet <- read.delim("gdc_sample_sheet.2023-06-27.tsv")
sample.sheet$Sample.ID <- gsub("-", "_", sample.sheet$Sample.ID)
n <- sample.sheet$Sample.ID

dir_res <- "../../../datos/ot/anakamura/scRNA_npj_v3/"

#Metadata of sample i
for(i in 1:length(n)){
  
print(paste0("Filtering", n[i]))
  
dir.create(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/FilteredByMorans"))
dir_filt <- paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/FilteredByMorans")
  
MDATA <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/metadata_", n[i], ".csv"), row.names = 1)

MDATA_AS <- MDATA[MDATA$singleR.labels.main == "Astrocyte",]
print(dim(MDATA_AS))
MDATA_MO <- MDATA[MDATA$singleR.labels.main %in% c("Macrophage", "Monocyte"),]
print(dim(MDATA_MO))
MDATA_T <- MDATA[MDATA$singleR.labels.main == "T_cells",]
print(dim(MDATA_T))

if(nrow(MDATA_T) != 0) { 
  #Load Morans stat
  
  DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/TCell_DEG_pr_test_res_", n[i], ".csv"))
  DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
  dim(DEG_morans)
  
  DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
  
  #Cluster j
  for(j in 1:length(unique(levels(as.factor(MDATA_T$seurat_clusters))))){
    cls <- as.numeric(unique(levels(as.factor(MDATA_T$seurat_clusters)))[j])
    
    #Expression matrix of cluster j
    
    EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_Counts_", n[i], ".tsv"),
                       row.names = 1)
    
    colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
    
    #MDATA of cluster j
    
    MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
    
    if(nrow(MDATA_j) == ncol(EXPM)){
      print(paste0("Cells check for cluster", cls) )
    } else { print(paste0("Cells disagree for cluster", cls)) }
    
    EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
    dim(EXPM)
    EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
    
    write.table(EXPM,
                file = paste0(dir_filt, "/TCells_Cls_", cls, "_MoransFilt", "_Counts_", n[i], ".tsv"),
                quote = FALSE,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE) 
  }  } else {
    print(paste0(n[i], " has no T cells"))
  }


if(nrow(MDATA_MO) != 0) { 
  #Load Morans stat
  
  DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/MO_DEG_pr_test_res_", n[i], ".csv"))
  DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
  dim(DEG_morans)
  
  DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
  
  #Cluster j
  for(j in 1:length(unique(levels(as.factor(MDATA_MO$seurat_clusters))))){
    cls <- as.numeric(unique(levels(as.factor(MDATA_MO$seurat_clusters)))[j])
    
    #Expression matrix of cluster j
    
    EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_Counts_", n[i], ".tsv"),
                       row.names = 1)
    
    colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
    
    #MDATA of cluster j
    
    MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
    
    if(nrow(MDATA_j) == ncol(EXPM)){
      print(paste0("Cells check for cluster", cls) )
    } else { print(paste0("Cells disagree for cluster", cls)) }
    
    EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
    dim(EXPM)
    EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
    
    write.table(EXPM,
                file = paste0(dir_filt ,"/MO_Cls_", cls, "_MoransFilt", "_Counts_", n[i], ".tsv"),
                quote = FALSE,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE) 
  } } else {
    print(paste0(n[i], " has no MO cells"))
  } 

if(nrow(MDATA_AS) != 0) { 
  #Load Morans stat
  
  DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/Astrocyte_DEG_pr_test_res_", n[i], ".csv"))
  DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
  dim(DEG_morans)
  
  DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
  
  #Cluster j
  for(j in 1:length(unique(levels(as.factor(MDATA_AS$seurat_clusters))))){
    cls <- as.numeric(unique(levels(as.factor(MDATA_AS$seurat_clusters)))[j])
    
    #Expression matrix of cluster j
    
    EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_Counts_", n[i], ".tsv"),
                       row.names = 1)
    
    colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
    
    #MDATA of cluster j
    
    MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
    
    if(nrow(MDATA_j) == ncol(EXPM)){
      print(paste0("Cells check for cluster", cls) )
    } else { print(paste0("Cells disagree for cluster", cls)) }
    
    EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
    dim(EXPM)
    EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
    
    write.table(EXPM,
                file = paste0(dir_filt, "/Astrocyte_Cls_", cls, "_MoransFilt", "_Counts_", n[i], ".tsv"),
                quote = FALSE,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE) 
  } } else {
    print(paste0(n[i], " has no AS cells"))
  } 
}

#For Scale.data
#Metadata of sample i
for(i in 1:length(n)){
  
  print(paste0("Filtering", n[i]))
  
  dir.create(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/FilteredByMorans"))
  dir_filt <- paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/FilteredByMorans")
  
  MDATA <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/metadata_", n[i], ".csv"), row.names = 1)
  
  MDATA_AS <- MDATA[MDATA$singleR.labels.main == "Astrocyte",]
  print(dim(MDATA_AS))
  MDATA_MO <- MDATA[MDATA$singleR.labels.main %in% c("Macrophage", "Monocyte"),]
  print(dim(MDATA_MO))
  MDATA_T <- MDATA[MDATA$singleR.labels.main == "T_cells",]
  print(dim(MDATA_T))
  
  if(nrow(MDATA_T) != 0) { 
    #Load Morans stat
    
    DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/TCell_DEG_pr_test_res_", n[i], ".csv"))
    DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
    dim(DEG_morans)
    
    DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
    
    #Cluster j
    for(j in 1:length(unique(levels(as.factor(MDATA_T$seurat_clusters))))){
      cls <- as.numeric(unique(levels(as.factor(MDATA_T$seurat_clusters)))[j])
      
      #Expression matrix of cluster j
      
      EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_ScaleData_", n[i], ".tsv"),
                         row.names = 1)
      
      colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
      
      #MDATA of cluster j
      
      MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
      
      if(nrow(MDATA_j) == ncol(EXPM)){
        print(paste0("Cells check for cluster", cls) )
      } else { print(paste0("Cells disagree for cluster", cls)) }
      
      EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
      dim(EXPM)
      EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
      
      write.table(EXPM,
                  file = paste0(dir_filt, "/TCells_Cls_", cls, "_MoransFilt", "_ScaleData_", n[i], ".tsv"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = TRUE,
                  row.names = FALSE) 
    }  } else {
      print(paste0(n[i], " has no T cells"))
    }
  
  
  if(nrow(MDATA_MO) != 0) { 
    #Load Morans stat
    
    DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/MO_DEG_pr_test_res_", n[i], ".csv"))
    DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
    dim(DEG_morans)
    
    DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
    
    #Cluster j
    for(j in 1:length(unique(levels(as.factor(MDATA_MO$seurat_clusters))))){
      cls <- as.numeric(unique(levels(as.factor(MDATA_MO$seurat_clusters)))[j])
      
      #Expression matrix of cluster j
      
      EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_ScaleData_", n[i], ".tsv"),
                         row.names = 1)
      
      colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
      
      #MDATA of cluster j
      
      MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
      
      if(nrow(MDATA_j) == ncol(EXPM)){
        print(paste0("Cells check for cluster", cls) )
      } else { print(paste0("Cells disagree for cluster", cls)) }
      
      EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
      dim(EXPM)
      EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
      
      write.table(EXPM,
                  file = paste0(dir_filt ,"/MO_Cls_", cls, "_MoransFilt", "_ScaleData_", n[i], ".tsv"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = TRUE,
                  row.names = FALSE) 
    } } else {
      print(paste0(n[i], " has no MO cells"))
    } 
  
  if(nrow(MDATA_AS) != 0) { 
    #Load Morans stat
    
    DEG_morans <- read.csv(paste0(dir_res, n[i],"/", n[i], "_results/TrajectoryAnalysis/Astrocyte_DEG_pr_test_res_", n[i], ".csv"))
    DEG_morans <- DEG_morans[DEG_morans$q_value < 1e-08,]
    dim(DEG_morans)
    
    DEG_morans <- DEG_morans %>% mutate(EnsemblNoV = annot_TCGA$Ensembl_noVers[match(DEG_morans$X, annot_TCGA$gene_id)])
    
    #Cluster j
    for(j in 1:length(unique(levels(as.factor(MDATA_AS$seurat_clusters))))){
      cls <- as.numeric(unique(levels(as.factor(MDATA_AS$seurat_clusters)))[j])
      
      #Expression matrix of cluster j
      
      EXPM <- read.delim(paste0(dir_res, n[i],"/", n[i], "_results/expression_matrices/Cls_", cls, "_ScaleData_", n[i], ".tsv"),
                         row.names = 1)
      
      colnames(EXPM) <- gsub("\\.", "-", colnames(EXPM))
      
      #MDATA of cluster j
      
      MDATA_j <- MDATA[rownames(MDATA) %in% colnames(EXPM),]
      
      if(nrow(MDATA_j) == ncol(EXPM)){
        print(paste0("Cells check for cluster", cls) )
      } else { print(paste0("Cells disagree for cluster", cls)) }
      
      EXPM <- EXPM[rownames(EXPM) %in% DEG_morans$EnsemblNoV,]
      dim(EXPM)
      EXPM <- EXPM %>% mutate(Gene = rownames(EXPM)) %>% dplyr::relocate(Gene)
      
      write.table(EXPM,
                  file = paste0(dir_filt, "/Astrocyte_Cls_", cls, "_MoransFilt", "_ScaleData_", n[i], ".tsv"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = TRUE,
                  row.names = FALSE) 
    } } else {
      print(paste0(n[i], " has no AS cells"))
    } 
}
