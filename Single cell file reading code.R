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
library(plyr)


################ Read the data and construct a Seurat object ################
data_dir <- paste0(getwd(), "/data") # Define the data storage path
samples = list.files(data_dir)       # List the sample files in the data folder
dir = file.path(data_dir, samples)   # Construct the file paths
afdata <- Read10X(data.dir = dir)

# Create a Seurat object
yjsl <- CreateSeuratObject(counts = afdata,
                           project = "SeuratObject",
                           min.cells = 3, # Genes must be expressed in at least 3 cells
                           min.features = 200) # At least 200 genes must be detected in each cell

# Assign cell identity information to the Type column for subsequent processing
afidens = mapvalues(Idents(yjsl), from = levels(Idents(yjsl)), to = samples)
Idents(yjsl) = afidens
yjsl$Type = Idents(yjsl)

################ Quality control visualization ################
# Calculate the percentage of mitochondrial genes
yjsl[["percent.mt"]] <- PercentageFeatureSet(yjsl, pattern = "^MT-")
# Calculate the percentage of ribosomal genes
yjsl[["percent.rb"]] <- PercentageFeatureSet(yjsl, pattern = "^RP")
# Draw violin plots showing the number of genes (nFeature_RNA), the number of detected molecules (nCount_RNA), and the percentage of mitochondrial genes (percent.mt) in each cell
VlnPlot(yjsl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3) #,pt.size = 0

## nCount represents the number of molecules detected in each cell
## nFeature represents the number of genes detected in each cell
### Outlier cells in either of the above metrics may be doublets or multiplets (the purpose of quality control is to remove outliers)
## Draw scatter plots to compare relationships between different features
plot2 <- FeatureScatter(yjsl, feature1 = "nCount_RNA", feature2 = "percent.rb") + RotatedAxis()
plot2
plot3 <- FeatureScatter(yjsl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis()
plot3

plot2 + plot3

### Filter genes
sub1 <- yjsl$nCount_RNA >= 1000 # Filter out cells with fewer than 1000 detected molecules
sub2 <- yjsl$nFeature_RNA >= 200 & yjsl$nFeature_RNA <= 10000 # Filter out cells with fewer than 200 or more than 10000 genes
sub3 <- yjsl$percent.mt <= 20 # Filter out cells with a mitochondrial gene percentage greater than 20%
sub4 <- yjsl$percent.rb <= 20 # Filter out cells with a ribosomal gene percentage greater than 20%
sub <- sub2 & sub3 & sub4 # Combine the filtering criteria
yjsl <- yjsl[, sub] # Apply the filtering criteria and retain cells that meet them

# Redraw the plots
VlnPlot(yjsl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

#### Official website
## Filter out cells with more than 2,500 or fewer than 200 genes
# Filter out cells with a mitochondrial gene percentage greater than 5%
yjsl <- subset(yjsl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Save for subsequent analysis
saveRDS(yjsl, "af.rds")