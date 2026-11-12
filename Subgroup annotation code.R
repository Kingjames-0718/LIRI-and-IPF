{
  library(dplyr)
  library(Seurat) # Single-cell data analysis
  library(patchwork) # Combine multiple ggplot plots into a single figure
  library(ggplot2)
  library(CCA)
  library(clustree)
  library(cowplot)
  library(monocle)
  library(tidyverse)
  library(SCpubr)
  library(GSEABase)
  library(harmony)
  library(plyr)
}

# Single-cell subpopulation annotation

load("yjsl_Type.rda")

## First, take a general look at the data
DimPlot(yjsl, reduction = "umap", label = T, repel = T)

## Check the specific cell groups
table(yjsl@meta.data[["cellType"]])

## Extract the subpopulation of interest
{
  Macrophage <- subset(yjsl, cellType == "Macrophage cells")
  # You can also extract a specific cluster using seurat_clusters == c("0", "9")
}

rm(yjsl)

## The subsequent steps follow the usual analysis workflow, including dimensionality reduction, clustering, and annotation
# Normalization
{
  Macrophage <- NormalizeData(
    Macrophage, normalization.method = "LogNormalize",
    scale.factor = 1e4
  )
  
  # Identify highly variable genes
  Macrophage <- FindVariableFeatures(
    Macrophage, selection.method = 'vst', nfeatures = 2000
  )
  
  # Scale the data as a necessary preparation for dimensionality reduction
  all.genes <- rownames(Macrophage)
  Macrophage <- ScaleData(Macrophage, features = all.genes)
  
  # Begin dimensionality reduction
  Macrophage <- RunPCA(
    Macrophage, features = VariableFeatures(object = Macrophage)
  )
}

# Data visualization methods
#VizDimReduction
VizDimLoadings(Macrophage, dims = 1:2, reduction = "pca")

#ElbowPlot
ElbowPlot(Macrophage, ndims = 20, reduction = "pca")

######## Determine the dimensions and cluster the cells ########
# Select PCs
PC = 1:15

### Cluster cells using the appropriate dimensions identified in the previous step
Macrophage <- FindNeighbors(Macrophage, dims = PC)

# resolution can be set between 0.1 and 1. Higher values produce more subclusters. A typical value is 0.5
Macrophage <- FindClusters(Macrophage, resolution = seq(0.2, 1.2, 0.1))
clustree(Macrophage)

Macrophage <- FindClusters(Macrophage, resolution = 0.5) ## Remember this dimension setting

### Perform nonlinear dimensionality reduction using UMAP
Macrophage <- RunUMAP(Macrophage, dims = PC)

# Data visualization
DimPlot(Macrophage, reduction = "umap", label = T, repel = T)

### Perform nonlinear dimensionality reduction using t-SNE
Macrophage = RunTSNE(Macrophage, dims = PC)
DimPlot(Macrophage, reduction = "tsne")

# Check how many types there are
table(Macrophage@meta.data$seurat_clusters)

# Retain only upregulated differentially expressed genes
markers <- FindAllMarkers(
  Macrophage, only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, file = "Macrophage_subset_marker_genes.csv")
head(markers)

# Annotation
{
  Macrophage <- FindClusters(Macrophage, resolution = 0.5)
  ## Remember to use the same resolution value as above
}

genes <- list("Macrophages" = c("CD163", "CD14"))

# Bubble plot
DotPlot(Macrophage, features = genes, cols = "RdYlBu") +
  RotatedAxis()

# Annotation
table(Macrophage@active.ident)
ann.ids <- c(
  "CD163+low", # cluster0
  "CD163+low", # cluster1
  "CD163+low", # Continue in order below
  "CD163+High",
  "CD163+low",
  "CD163+low",
  "CD163+low",
  "CD163+low"
)

afidens = mapvalues(
  Idents(Macrophage),
  from = levels(Idents(Macrophage)),
  to = ann.ids
)
Idents(Macrophage) = afidens
Macrophage$newType = Idents(Macrophage)

## Take another look

DimPlot(Macrophage, reduction = "umap", label = T, label.size = 3.5) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black", size = 0.5, linetype = "solid"
    ),
    legend.position = "right"
  )

library(ggsci)
cors <- pal_igv()(16) # Define colors

DimPlot(
  Macrophage, reduction = "umap", label = T,
  group.by = "newType", cols = cors
)

# Save the file with assigned cell types for subsequent cell-cell communication and pseudotime trajectory analyses
save(Macrophage, file = "Macrophage.rda")