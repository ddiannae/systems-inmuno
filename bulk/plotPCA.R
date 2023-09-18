library(ggplot2)
library(tidyr)

# Raw counts PCA, plotPCA code modified from DESeq2::plotPCA 
# getMethod("plotPCA", "DESeqTransform")
plotPCA = function(counts, data, intgroup="group", ntop = 500, pcsToUse=1:2)
{
  # calculate the variance for each gene
  rv <- rowVars(counts, useNames = TRUE)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(counts[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2 )
  
  if (!all(intgroup %in% names(data))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(data[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    data[[intgroup]]
  }
  
  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, intgroup.df, name=colnames(counts))
  colnames(d)[1:2] <- pcs
  
ggplot(data=d, aes(x=.data[[pcs[1]]], y=.data[[pcs[2]]], color=.data[[intgroup]])) +
    geom_point(size=3) + 
    xlab(paste0(pcs[1],": ",round(percentVar[pcsToUse[1]] * 100),"% variance")) +
    ylab(paste0(pcs[2],": ",round(percentVar[pcsToUse[2]] * 100),"% variance")) +
    coord_fixed()  + 
    scale_color_brewer(palette = "Dark2") +
    theme_bw()
}


plotBoxplot = function(counts, data, intgroup="group", id = "id", type = "sample", scaled = TRUE)
{
  if(scaled) {
    d <- pivot_longer(as_tibble(t(scale(t(counts), center = FALSE)), rownames = "gene"), cols = !matches("gene"),
                      names_to = id)
    yname <- "Scaled expression"
  } else {
    d <- pivot_longer(as_tibble(counts, rownames = "gene"), cols = !matches("gene"),
                      names_to = id)
    yname <- "Expression"
  }
  
  d <- d %>% inner_join(data, by = id)
  
  if(type == "sample") {
    ggplot(data=d, aes(y = value, x =.data[[id]],  color=.data[[intgroup]])) +
      geom_boxplot(notch = TRUE) + 
      ylab(yname) + 
      scale_color_brewer(palette = "Dark2") +
      theme_classic()
  } else if (type == "group") {
    ggplot(data=d, aes(y = value, x =.data[[intgroup]],  color=.data[[intgroup]])) +
      geom_boxplot(notch = TRUE) +
      ylab(yname) +
      scale_color_brewer(palette = "Dark2") +
      theme_classic()
  }
}