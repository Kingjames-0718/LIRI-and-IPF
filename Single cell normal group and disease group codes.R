library(dplyr)
library(Seurat) # Single-cell data analysis
library(patchwork) # Combine multiple ggplot plots into a single figure
library(ggplot2)
library(SingleR) # Automatically annotate single-cell data
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(GSEABase)
library(harmony)


################ Read the data and construct a Seurat object ################
{
  data_dir <- paste0(getwd(), "/data") # Define the data storage path
  samples = list.files(data_dir)      # List the sample files in the data folder
  dir = file.path(data_dir, samples)  # Construct the file paths
  afdata <- Read10X(data.dir = dir)
  
  # Create a Seurat object
  yjsl <- CreateSeuratObject(
    counts = afdata,
    project = "SeuratObject",
    min.cells = 3, # Genes must be expressed in at least 3 cells
    min.features = 200 # At least 200 genes must be detected in each cell
  )
  
  # Assign cell identity information to the Type column for subsequent processing
  afidens = mapvalues(Idents(yjsl), from = levels(Idents(yjsl)), to = samples)
  Idents(yjsl) = afidens
  yjsl$Type = Idents(yjsl)
}


# Define a new group variable based on the sample names in Type
yjsl$group <- ifelse(
  yjsl$Type %in% c("GSM4506698", "GSM4506700"), "NC",
  ifelse(yjsl$Type %in% c("GSM4506699", "GSM4506701"), "Tumor", NA)
)

# View the first few rows to confirm the grouping information
head(yjsl@meta.data)


{
  ################ Quality control visualization ################
  # Calculate the percentage of mitochondrial genes
  yjsl[["percent.mt"]] <- PercentageFeatureSet(yjsl, pattern = "^MT-")
  
  # Calculate the percentage of ribosomal genes
  yjsl[["percent.rb"]] <- PercentageFeatureSet(yjsl, pattern = "^RP")
  
  # Draw violin plots showing the number of genes (nFeature_RNA), the number of detected molecules (nCount_RNA), and the percentage of mitochondrial genes (percent.mt) in each cell
  VlnPlot(
    yjsl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  ) #,pt.size = 0
  
  ## nCount represents the number of molecules detected in each cell
  ## nFeature represents the number of genes detected in each cell
  ### Outlier cells in either of these metrics may be doublets or multiplets (the purpose of quality control is to remove outliers)
  
  ## Draw scatter plots to compare relationships between different features
  plot2 <- FeatureScatter(
    yjsl, feature1 = "nCount_RNA", feature2 = "percent.rb"
  ) + RotatedAxis()
  plot2
  
  plot3 <- FeatureScatter(
    yjsl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
  ) + RotatedAxis()
  plot3
  
  plot2 + plot3
  
  ### Filter genes
  sub1 <- yjsl$nCount_RNA >= 1000 # Filter out cells with fewer than 1000 detected molecules
  sub2 <- yjsl$nFeature_RNA >= 200 & yjsl$nFeature_RNA <= 10000 # Filter out cells with fewer than 200 or more than 10000 genes
  sub3 <- yjsl$percent.mt <= 20 # Filter out cells with a mitochondrial gene percentage greater than 20%
  sub4 <- yjsl$percent.rb <= 20 # Filter out cells with a ribosomal gene percentage greater than 20%
  sub <- sub2 & sub3 & sub4 # Combine the filtering criteria
  yjsl <- yjsl[, sub] # Apply the filtering criteria and retain cells that meet them
}

# Redraw the plots
VlnPlot(
  yjsl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)


####### Use the LogNormalize method #######
yjsl <- NormalizeData(
  yjsl, normalization.method = "LogNormalize",
  scale.factor = 10000
)

####### Identify highly variable genes #######

yjsl <- FindVariableFeatures(yjsl, selection.method = "vst", nfeatures = 2000)

# Extract the top 15 highly variable genes
top15 <- head(VariableFeatures(yjsl), 15)
top15

## Display highly variable genes
plot1 <- VariableFeaturePlot(yjsl, pt.size = 2, raster = T)
plot1

plot2 <- LabelPoints(
  plot = plot1, points = top15,
  xnudge = 0,
  ynudge = 0
)
plot2


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

plot2 <- LabelPoints(
  plot = plot1, points = top15,
  xnudge = 0,
  ynudge = 0
)
plot2


###### Scale the data as a necessary preparation for dimensionality reduction
all.genes <- rownames(yjsl)
yjsl <- ScaleData(yjsl, features = all.genes)

#### Begin dimensionality reduction ###
yjsl <- RunPCA(yjsl, features = VariableFeatures(object = yjsl))

# Data visualization methods

#DimHeatmap
DimHeatmap(yjsl, dims = 1:15, cells = 500, balanced = TRUE)

# ElbowPlot (the elbow point can be used to determine the number of dimensions)
ElbowPlot(yjsl, ndims = 20, reduction = "pca")


######## Determine the dimensions and cluster the cells ########
# Select PCs
PC = 1:15

### Cluster cells using the appropriate dimensions identified in the previous step
yjsl <- FindNeighbors(yjsl, dims = PC)

# resolution can be set between 0.1 and 1. Higher values produce more subclusters. A typical value is 0.5
yjsl <- FindClusters(yjsl, resolution = seq(0.2, 1.2, 0.1))
clustree(yjsl)


## !!! Remember this dimension setting

yjsl <- FindClusters(yjsl, resolution = 0.5)

### Perform nonlinear dimensionality reduction using UMAP
yjsl <- RunUMAP(yjsl, dims = PC)

# Data visualization
DimPlot(yjsl, reduction = "umap", label = T, repel = T)

DimPlot(yjsl, reduction = "umap", label = T, repel = T, split.by = 'group')

### Perform nonlinear dimensionality reduction using t-SNE
yjsl = RunTSNE(yjsl, dims = PC)
embed_tsne <- Embeddings(yjsl, 'tsne')
DimPlot(yjsl, reduction = "tsne")
DimPlot(yjsl, reduction = "tsne", split.by = 'group')


### Perform differential expression analysis between the normal and tumor groups
# Assign the group column as the identities of the Seurat object
## Remember to run this line frequently during subsequent analyses
Idents(yjsl) <- yjsl$group

# Perform differential expression analysis
{
  DEGs <- FindMarkers(
    yjsl, ident.1 = 'Tumor', ident.2 = 'NC',
    logfc.threshold = 0.58
  ) # Adjustable
}

# Add grouping information to the differentially expressed gene table
DEGs$group <- ifelse(DEGs$avg_log2FC > 0, "Tumor", "NC") # Determine the group based on log2FoldChange

# Export the differential expression results
write.csv(DEGs, "differentially_expressed_genes_tumor_vs_normal.csv")

## Filter by p-values and adjusted p-values
DEGs.filter <- subset(DEGs, p_val_adj < 0.05)
write.csv(DEGs.filter, "filtered_Tumor_vs_NC_differentially_expressed_genes.csv")

## Restore the original identities
Idents(yjsl) <- yjsl$seurat_clusters

### Take a general look at the results ###
# Check how many types there are
table(yjsl@meta.data$seurat_clusters)

# Retain only upregulated differentially expressed genes
markers <- FindAllMarkers(
  yjsl, only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, file = "cluster_marker_genes.csv")
head(markers)

## Save for subsequent analysis
# saveRDS(yjsl, "af2.rds")

library(randomcoloR)

## A package specifically for enhancing plot appearance
#install.packages("/Users/zhangzhenhu/Desktop/single_cell_annotation_manual/scRNAtoolVis_0.1.0.tar.gz", repos=NULL, type="source")
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

# yjsl = readRDS("af2.rds")

################ Recommended method of manual annotation ################
###### First approach using differentially expressed genes ################

# Identify all differentially expressed genes
Idents(yjsl) <- yjsl$seurat_clusters
allmarkers <- FindAllMarkers(
  yjsl, only.pos = TRUE,
  min.pct = 0.25, # Only genes expressed in at least 25% of cells are considered for differential expression
  logfc.threshold = 0.25 # Select only genes with a log2 fold change greater than 0.25
)


# save(allmarkers, file = "allmarkers.rda")
# load("allmarkers.rda")

# Extract the top 4 genes with the most significant differential expression in each cluster
top4_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC)


# Heatmap with a different color scheme
DoHeatmap(yjsl, features = unique(top4_markers$gene)) + NoLegend() +
  scale_fill_gradientn(colors = c("#2fa1dd", "white", "#f87669"))

# Bubble plot
DotPlot(yjsl, features = unique(top4_markers$gene), cols = "RdYlBu") +
  RotatedAxis()

### You can then visit annotation websites
# CellMarker http://xteam.xbio.top/CellMarker/

###### Second approach using common marker genes ################
# T Cells (CD3D, CD3E, CD8A)
#
# B cells (CD19, CD79A, MS4A1 [CD20])
#
# Plasma cells (IGHG1, MZB1, SDC1, CD79A)
#
# Monocytes and macrophages (CD68, CD163, CD14)
#
# NK Cells (FGFBP2, FCG3RA, CX3CR1)
#
# Photoreceptor cells (RCVRN)
#
# Fibroblasts (FGF7, MME)
#
# Endothelial cells (PECAM1, VWF)
#
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24)
#
# immune (CD45+,PTPRC)
#
# epithelial/cancer (EpCAM+,EPCAM)
#
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo)
#
# mast cells( TPSAB1 and TPSB2 )
#
# naive B cells(MS4A1 (CD20), CD19, CD22, TCL1A, and CD83)
#
# plasma B cells(CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4


####### Consult published studies ################

yjsl <- FindClusters(yjsl, resolution = 0.5) ## Remember to use the same resolution value as above

genes <- list(
  "T cells" = c("CD3D", "CD2"),
  "Macrophages" = c("CD163", "CD14"),
  "B cells" = c("CD79A", "MS4A1"),
  "Fibroblasts" = c("ACTA2", "PDGFRB"),
  "Endothelial cells" = c("VWF", "RAMP2"),
  "Mast cells" = c("TPSAB1", "TPSB2"),
  "Epithelial cells" = c("SFTPA1", "EPCAM")
)

do_DotPlot(
  sample = yjsl, features = genes, dot.scale = 10,
  colors.use = c("yellow", "red"), legend.length = 50,
  legend.framewidth = 2, font.size = 10
)

# Bubble plot with another color scheme
DotPlot(yjsl, features = genes, cols = "RdYlBu") +
  RotatedAxis()

# Manual annotation
table(yjsl@active.ident)
ann.ids <- c(
  "Epithelial cells", # cluster0
  "T cells",          # cluster1
  "Fibroblasts",      # Continue in order below
  "T cells",
  "Macrophage cells",
  "Endothelial cells",
  "Macrophage cells",
  "Epithelial cells",
  "Epithelial cells",
  "Macrophage cells",
  "Mast cells",
  "Epithelial cells",
  "B cells",
  "Mast cells",
  "Epithelial cells",
  "Fibroblasts"
)

afidens = mapvalues(Idents(yjsl), from = levels(Idents(yjsl)), to = ann.ids)
Idents(yjsl) = afidens
yjsl$cellType = Idents(yjsl)

######### Visualize the results after manual annotation
# Visualize UMAP/t-SNE

DimPlot(yjsl, reduction = "umap", label = T, label.size = 3.5) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black", size = 0.5, linetype = "solid"
    ),
    legend.position = "right"
  )

DimPlot(yjsl, reduction = "tsne", label = T, label.size = 3.5) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black", size = 0.5, linetype = "solid"
    ),
    legend.position = "right"
  )