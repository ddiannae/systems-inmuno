library(monocle3)
library(monocle)
library(Seurat)
library(SeuratWrappers)
library(SummarizedExperiment)
library(dplyr)

#Load annotation
annot_TCGA <- readRDS("row_data.RDS")

#Load Seurat object
files <- list.files(pattern = "seurat.*RDS")

#Get markers and cluster
#Get sample name
n <- gsub("^seurat_", "", files)
n <- gsub("\\.RDS$", "", n)
n

##### Get outs ####
for(i in 1:length(n)){
  
  MD <- read.delim(paste0("../../../datos/ot/anakamura/scRNA_npj/", n[i], "_results", "/metadata_", n[i],".csv"),
                   sep = ",")
  
  markers <- read.delim(paste0("../../../datos/ot/anakamura/scRNA_npj/", n[i], "_results", "/markers_", n[i],".tsv"),
                        sep = " ")
  
  markers$HGNC <- annot_TCGA$gene_name[match(markers$gene, annot_TCGA$gene_id)]
  
  marker_genes_AS <- c("CD44", "CD133", "SOX2", "KLF4", "POU5F1", "CD90", "CD15",
                    "CXCR4", "CX3CR1", "MYC", "NANOG", "NES", "MSI1")
  
  marker_genes_MO <- c("CTSB", "CD163")
  
  marker_genes_T <- c("B2M", "INPP4B")

  clusters_AS <- MD[MD$singleR.labels.main == "Astrocyte",]$seurat_clusters
  
  if(length(clusters_AS) == 0) {
    print(paste0("Sample", n[i], "has no AS cells"))
  } else {
    
    markers_sub_AS <- markers[markers$cluster %in% unique(clusters_AS),]
    
    if(exists("outs_AS") == FALSE){
      outs_AS <- list(markers_sub_AS[markers_sub_AS$HGNC %in% marker_genes_AS,])
    } else {
      AS <- markers_sub_AS[markers_sub_AS$HGNC %in% marker_genes_AS,]
      outs_AS <- c(outs_AS, list(AS))
    }

  }
    
  clusters_MO <- MD[MD$singleR.labels.main %in% c("Macrophage", "Monocyte"),]$seurat_clusters
  
  if(length(clusters_MO) == 0) {
    print(paste0("Sample", n[i], "has no MO cells"))
  } else {
    
    markers_sub_MO <- markers[markers$cluster %in% clusters_MO,]
    
    if(exists("outs_MO") == FALSE){
      outs_MO <- list(markers_sub_MO[markers_sub_MO$HGNC %in% marker_genes_MO,])
    } else {
      MO <- markers_sub_MO[markers_sub_MO$HGNC %in% marker_genes_MO,]
      outs_MO <- c(outs_MO, list(MO))
    }
    
  }
  
  clusters_T <- MD[MD$singleR.labels.main == "T_cells",]$seurat_clusters
  
  if(length(clusters_T) == 0) {
    print(paste0("Sample ", n[i], " has no T cells"))
  } else {
    
    markers_sub_T <- markers[markers$cluster %in% clusters_T,]
    
    if(exists("outs_T") == FALSE){
      outs_T <- list(markers_sub_T[markers_sub_T$HGNC %in% marker_genes_T,])
    } else {
      Tc <- markers_sub_T[markers_sub_T$HGNC %in% marker_genes_T,]
      outs_T <- c(outs_T, list(Tc))
    }
  }
}

# [1] "SampleCPT0087680014has no T cells"
# [1] "Sample CPT0125220004 has no T cells"
# [1] "SampleCPT0167860015has no T cells"
# [1] "SampleCPT0189650015has no T cells"
# [1] "SampleCPT0207030018has no T cells"
# [1] "SampleCPT0224600013has no T cells"  

NoTcelss <- c("CPT0087680014", "CPT0125220004", "CPT0167860015",
              "CPT0189650015", "CPT0207030018", "CPT0224600013")

  names(outs_AS) <- n
  names(outs_MO) <- n
  names(outs_T) <- n[!(n %in% NoTcelss)]
  
#### Get TRajectories #### 
  
for(i in 1:length(n)){
  
  
  dir.create(paste0("../../../datos/ot/anakamura/scRNA_npj/",n[i], "_results", "/TrajectoryAnalysis"))
  dir_res <- paste0("../../../datos/ot/anakamura/scRNA_npj/",n[i], "_results", "/TrajectoryAnalysis")
  
  #Load seurat
  
  seurat.obj <- readRDS(paste0("../../../datos/ot/anakamura/scRNA_npj/",n[i], "_results", "/Processed_Seurat_", n[i], ".RDS"))
  
  ##### Subset Astrocytes ####
  
  print(paste0("Analyzing ", n[i], " Astrocytes"))
  
  cells <- seurat.obj[,seurat.obj$singleR.labels.main == "Astrocyte"]
  
  #Transform to Cell Data Set
  cells.cds <- as.cell_data_set(cells)
  fData(cells.cds)
  fData(cells.cds)$gene_short_name <- rownames(fData(cells.cds))
  
  #Construct trajectory
  cells.cds <- cluster_cells(cds = cells.cds, reduction_method = "UMAP")
  cells.cds <- learn_graph(cells.cds, use_partition = TRUE)
  
  #Get root using the expression of stemness, priorizing SOX2 expression
  get_earliest_principal_node <- function(cds){
    
    table <- which(names(outs_AS) == n[i])
    
    if("SOX2" %in% outs_AS[[table]]$HGNC){
      
      if(length(c(outs_AS[[table]][outs_AS[[table]]$HGNC == "SOX2",]$cluster)) > 1){
        cls_no <- which(outs_AS[[table]][outs_AS[[table]]$HGNC == "SOX2",]$p_val_adj == min(outs_AS[[table]][outs_AS[[table]]$HGNC == "SOX2",]$p_val_adj))
        
        cls <- outs_AS[[table]][outs_AS[[table]]$HGNC == "SOX2",]$cluster[cls_no]
        
      } else {
        
        cls <- outs_AS[[table]][outs_AS[[table]]$HGNC == "SOX2",]$cluster
      }
      
      cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
      
      closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
      
    } else {
      
      if("CD44" %in% outs_AS[[table]]$HGNC) {
        
        cls_no <- which(outs_AS[[table]][outs_AS[[table]]$HGNC == "CD44",]$p_val_adj == min(outs_AS[[table]][outs_AS[[table]]$HGNC == "CD44",]$p_val_adj))
        
        cls <- outs_AS[[table]][outs_AS[[table]]$HGNC == "CD44",]$cluster[cls_no]
        
        cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
        
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
        
      } else {
        
        cls <-outs_AS[[table]][outs_AS[[table]]$HGNC == "NES" ,]$cluster
        
        cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
        
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
        
      }
    }
  }
  
  cells.cds <- order_cells(cells.cds, root_pr_nodes=get_earliest_principal_node(cells.cds))
  
  # plot trajectories colored by pseudotime
  trajectory <- plot_cells(
    cds = cells.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  
  ggsave(filename =  paste0(dir_res, "/Astrocyte_Trayectory_UMAP_", n[i], ".tiff"), plot = trajectory, units="in", 
         width=10, height=7,  dpi =300)

  #Get DEG
  cells_cds_pr_test_res <- graph_test(cells.cds, neighbor_graph="principal_graph", cores=70)
  cells_cds_pr_test_res$gene_short_name <- annot_TCGA$gene_name[match(rownames(cells_cds_pr_test_res), annot_TCGA$gene_id)]
  
  write.csv(cells_cds_pr_test_res, file = paste0(dir_res, "/Astrocyte_DEG_pr_test_res_", n[i], ".csv"), 
            quote = FALSE)
  
  #Plot DEG genes
 pr_deg_genes <- plot_cells(cells.cds, genes=rownames(cells_cds_pr_test_res %>% arrange(desc(morans_test_statistic)))[1:42],
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
  
  ggsave(filename =  paste0(dir_res, "/Astrocyte_pr_DEG_UMAP_", n[i], ".tiff"), plot = pr_deg_genes, units="in", 
         width=20, height=17,  dpi =300)
  
  #Save Cell Data Set obj
  saveRDS(cells.cds, paste0(dir_res, "/Cell_Data_Set_Astrocytes",n[i],".RDS"))
  
  ##### Subset MO ####
  
  print(paste0("Analyzing ", n[i], " MO"))
  
  cells <- seurat.obj[,seurat.obj$singleR.labels.main %in% c("Macrophage", "Monocyte")]
  
  #Transform to Cell Data Set
  cells.cds <- as.cell_data_set(cells)
  fData(cells.cds)
  fData(cells.cds)$gene_short_name <- rownames(fData(cells.cds))
  
  #Construct trajectory
  cells.cds <- cluster_cells(cds = cells.cds, reduction_method = "UMAP")
  cells.cds <- learn_graph(cells.cds, use_partition = TRUE)
  
  #Get root using the expression of stemness, priorizing SOX2 expression
  get_earliest_principal_node <- function(cds){
    
    if("Monocyte" %in% seurat.obj$singleR.labels.main) {
      
      cell_ids <- which(colData(cells.cds)[, "singleR.labels.main"] == "Monocyte")
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    } else {
      
      table <- which(names(outs_MO) == n[i])
      
      cls_no <- which(outs_MO[[table]][outs_MO[[table]]$HGNC == "CD163",]$p_val_adj == min(outs_MO[[table]][outs_MO[[table]]$HGNC == "CD163",]$p_val_adj))
      
      cls <- outs_MO[[table]][outs_MO[[table]]$HGNC == "CD163",]$cluster[cls_no]
      
      cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
      
    }
  }
  
  cells.cds <- order_cells(cells.cds, root_pr_nodes=get_earliest_principal_node(cells.cds))
  
  # plot trajectories colored by pseudotime
  trajectory <- plot_cells(
    cds = cells.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )

  cells.cds <- order_cells(cells.cds, root_pr_nodes=get_earliest_principal_node(cells.cds))
  
  # plot trajectories colored by pseudotime
  trajectory <- plot_cells(
    cds = cells.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  
  
  ggsave(filename =  paste0(dir_res, "/MO_Trayectory_UMAP_", n[i], ".tiff"), plot = trajectory, units="in", 
         width=10, height=7,  dpi =300)
  
  #Get DEG
  cells_cds_pr_test_res <- graph_test(cells.cds, neighbor_graph="principal_graph", cores=70)
  cells_cds_pr_test_res$gene_short_name <- annot_TCGA$gene_name[match(rownames(cells_cds_pr_test_res), annot_TCGA$gene_id)]
  
  write.csv(cells_cds_pr_test_res, file = paste0(dir_res, "/MO_DEG_pr_test_res_", n[i], ".csv"), 
            quote = FALSE)
  
  #Plot DEG genes
  pr_deg_genes <- plot_cells(cells.cds, genes=rownames(cells_cds_pr_test_res %>% arrange(desc(morans_test_statistic)))[1:42],
                             show_trajectory_graph=FALSE,
                             label_cell_groups=FALSE,
                             label_leaves=FALSE)
  
  ggsave(filename =  paste0(dir_res, "/MO_pr_DEG_UMAP_", n[i], ".tiff"), plot = pr_deg_genes, units="in", 
         width=20, height=17,  dpi =300)
  
  #Save Cell Data Set obj
  saveRDS(cells.cds, paste0(dir_res, "/Cell_Data_Set_MO_",n[i],".RDS"))
  
  #### Subset T cells ####
  if(n[i] %in% NoTcelss) {
    print(paste0("Sample ", n[i], "has no T cells"))
  } else {
    
    print(paste0("Analyzing ", n[i], " T cells"))
    
    cells <- seurat.obj[,seurat.obj$singleR.labels.main == "T_cells"]
    
    #Transform to Cell Data Set
    cells.cds <- as.cell_data_set(cells)
    fData(cells.cds)
    fData(cells.cds)$gene_short_name <- rownames(fData(cells.cds))
    
    #Construct trajectory
    cells.cds <- cluster_cells(cds = cells.cds, reduction_method = "UMAP")
    cells.cds <- learn_graph(cells.cds, use_partition = TRUE)
    
    #Get root using the expression of stemness, priorizing SOX2 expression
    get_earliest_principal_node <- function(cds){
      
      table <- which(names(outs_T) == n[i])

      if(length(c(outs_T[[table]][outs_T[[table]]$HGNC == "B2M",]$cluster)) > 2) {
        
        cls_no <- which(outs_T[[table]][outs_T[[table]]$HGNC == "B2M",]$p_val_adj == min(outs_T[[table]][outs_T[[table]]$HGNC == "B2M",]$p_val_adj))
        
        cls <- outs_T[[table]][outs_T[[table]]$HGNC == "B2M",]$cluster[cls_no]
        
        cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
        
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
 
      } else {
        
        cls <- outs_T[[table]][outs_T[[table]]$HGNC == "B2M" ,]$cluster
        
        cell_ids <- which(colData(cds)[, "seurat_clusters"] == cls)
        
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
        
        
      }
    }
      
    cells.cds <- order_cells(cells.cds, root_pr_nodes=get_earliest_principal_node(cells.cds))
    
    # plot trajectories colored by pseudotime
    trajectory <- plot_cells(
      cds = cells.cds,
      color_cells_by = "pseudotime",
      show_trajectory_graph = TRUE
    )
    
    ggsave(filename =  paste0(dir_res, "/TCell_Trayectory_UMAP_", n[i], ".tiff"), plot = trajectory, units="in", 
           width=10, height=7,  dpi =300)
    
    #Get DEG
    cells_cds_pr_test_res <- graph_test(cells.cds, neighbor_graph="principal_graph", cores=70)
    cells_cds_pr_test_res$gene_short_name <- annot_TCGA$gene_name[match(rownames(cells_cds_pr_test_res), annot_TCGA$gene_id)]
    
    write.csv(cells_cds_pr_test_res, file = paste0(dir_res, "/TCell_DEG_pr_test_res_", n[i], ".csv"), 
              quote = FALSE)
    
    #Plot DEG genes
    pr_deg_genes <- plot_cells(cells.cds, genes=rownames(cells_cds_pr_test_res %>% arrange(desc(morans_test_statistic)))[1:42],
                               show_trajectory_graph=FALSE,
                               label_cell_groups=FALSE,
                               label_leaves=FALSE)
    
    ggsave(filename =  paste0(dir_res, "/TCell_pr_DEG_UMAP_", n[i], ".tiff"), plot = pr_deg_genes, units="in", 
           width=20, height=17,  dpi =300)
    
    
    #Save Cell Data Set obj
    saveRDS(cells.cds, paste0(dir_res, "/Cell_Data_Set_TCell_",n[i],".RDS"))
    
    
  }
  
}
  
