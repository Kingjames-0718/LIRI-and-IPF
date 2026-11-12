library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
#BiocManager::install("UCell")
library(UCell)
library(AUCell)
library(decoupleR)
library(singscore)
#devtools::install_github("chuiqin/irGSEA", force = T)
library(irGSEA) ## Installation is relatively complicated
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)

### Begin the analysis ###
yjsl = readRDS("af.rds")

####### Use the LogNormalize method #######
yjsl <- NormalizeData(yjsl, normalization.method = "LogNormalize",
                      scale.factor = 10000)

####### Identify highly variable genes #######

yjsl <- FindVariableFeatures(yjsl, selection.method = "vst", nfeatures = 2000)

# Extract the top 15 highly variable genes
top15 <- head(VariableFeatures(yjsl), 15)
top15

## Display highly variable genes
plot1 <- VariableFeaturePlot(yjsl, pt.size = 2, raster = T)
plot1
plot2 <- LabelPoints(plot = plot1, points = top15,
                     xnudge = 0,
                     ynudge = 0)
plot2
plot1 + plot2

###### Scale the data as a necessary preparation for dimensionality reduction
{
  all.genes <- rownames(yjsl)
  yjsl <- ScaleData(yjsl, features = all.genes)
}

#### Begin dimensionality reduction ###
yjsl <- RunPCA(yjsl, features = VariableFeatures(object = yjsl))
# Data visualization methods
#VizDimReduction
VizDimLoadings(yjsl, dims = 1:2, reduction = "pca")
# DimPlot
DimPlot(yjsl, reduction = "pca")
#DimHeatmap
DimHeatmap(yjsl, dims = 1:15, cells = 500, balanced = TRUE)
#ElbowPlot
ElbowPlot(yjsl, ndims = 20, reduction = "pca")

######## Remove batch effects using Harmony ########################
yjsl <- RunHarmony(yjsl, group.by.vars = "Type")

######## The following steps are the same ########

####### Identify highly variable genes #######

yjsl <- FindVariableFeatures(yjsl, selection.method = "vst", nfeatures = 2000)

# Extract the top 15 highly variable genes
top15 <- head(VariableFeatures(yjsl), 15)
top15

## Display highly variable genes
plot1 <- VariableFeaturePlot(yjsl, pt.size = 2, raster = T)
plot1
plot2 <- LabelPoints(plot = plot1, points = top15,
                     xnudge = 0,
                     ynudge = 0)
plot2
plot1 + plot2

###### Scale the data as a necessary preparation for dimensionality reduction
all.genes <- rownames(yjsl)
yjsl <- ScaleData(yjsl, features = all.genes)

#### Begin dimensionality reduction ###
yjsl <- RunPCA(yjsl, features = VariableFeatures(object = yjsl))
# Data visualization methods
#VizDimReduction
VizDimLoadings(yjsl, dims = 1:2, reduction = "pca")
# DimPlot
DimPlot(yjsl, reduction = "pca")
#DimHeatmap
DimHeatmap(yjsl, dims = 1:15, cells = 500, balanced = TRUE)

# ElbowPlot (the elbow point can be used to determine the number of dimensions)
ElbowPlot(yjsl, ndims = 20, reduction = "pca")


# ###### You can directly use the dimensions indicated by the scatter plot above. Either approach is acceptable, depending on your needs (official website) ####
# yjsl <- JackStraw(yjsl, num.replicate = 100) # Takes a long time
# yjsl <- ScoreJackStraw(yjsl, dims = 1:20)
# # Select the dimensions with the most significant changes in p-values
# JackStrawPlot(yjsl, dims = 1:15)
# # Determine the elbow value, which is 10 ± 1 in this example
# ElbowPlot(yjsl)

######## Determine the dimensions and cluster the cells ########
# Select PCs
PC = 1:10
### Cluster cells using the appropriate dimensions identified in the previous step
yjsl <- FindNeighbors(yjsl, dims = PC)
# resolution can be set between 0.1 and 1. Higher values produce more subclusters. A typical value is 0.5
yjsl <- FindClusters(yjsl, resolution = seq(0.2, 1.2, 0.1))
clustree(yjsl)

yjsl <- FindClusters(yjsl, resolution = 0.5) ## Remember this dimension setting
### Perform nonlinear dimensionality reduction using UMAP
yjsl <- RunUMAP(yjsl, dims = PC)
# Data visualization
DimPlot(yjsl, reduction = "umap", label = T, repel = T)


### Perform nonlinear dimensionality reduction using t-SNE
yjsl = RunTSNE(yjsl, dims = PC)
embed_tsne <- Embeddings(yjsl, 'tsne')
DimPlot(yjsl, reduction = "tsne")


### Take a general look at the results ###
# Check how many types there are
table(yjsl@meta.data$seurat_clusters)
# Retain only upregulated differentially expressed genes
markers <- FindAllMarkers(yjsl, only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers, file = "cluster_marker_genes.csv")
head(markers)

## Save for subsequent analysis
saveRDS(yjsl, "af2.rds")