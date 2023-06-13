library(Seurat)
library(SeuratDisk)
library(SummarizedExperiment)
library(TCGAbiolinks)

#Get barcodes

qry <- GDCquery(project = "CPTAC-3",
                data.category= "Transcriptome Profiling",
                data.type = "Single Cell Analysis",
                workflow.type = "Seurat - 10x Chromium",
                access = "open",
                data.format = "HDF5")

clinical <- read.delim("clinical.tsv")
qry.t <- qry[[1]][[1]]
qry.t <- qry.t[qry.t$cases.submitter_id %in% clinical$case_submitter_id,]

barcodes <- qry.t$cases
saveRDS(barcodes, "barcodes.RDS")
barcodes <- readRDS("barcodes.RDS")


#Get HDF5s

for(i in 1:length(barcodes)) {
  
  qry <- GDCquery(project = "CPTAC-3",
                  data.category= "Transcriptome Profiling",
                  data.type = "Single Cell Analysis",
                  workflow.type = "Seurat - 10x Chromium",
                  access = "open",
                  barcode = barcodes[i],
                  data.format = "HDF5")
  
  GDCdownload(qry)
  
  saveRDS(LoadLoom(paste0("GDCdata/CPTAC-3/harmonized/Transcriptome_Profiling/Single_Cell_Analysis/", 
                          qry.t[i,1],
                          "/seurat.loom")),
          paste0("seurat_", barcodes[i], ".RDS"))
  
}

