setwd("")

library(Seurat)
library(tidyverse)
library(ggplot2)
######
load("yjsl_Type.rda")

DimPlot(yjsl, group.by = 'cellType', label = TRUE)

# The key is to input your own dataset
# Custom gene set
Gene <- data.table::fread("gene_set.txt")

{
  gene <- Gene$gene
  ## The required format is a list
  genes = list(gene)
}

colnames(yjsl@meta.data)
# Core code for calculating the score
yjsl = AddModuleScore(yjsl, features = genes, name = 'autophagy')

## Take a general look at the distribution
plot(yjsl$autophagy1)

# Visualize the results
FeaturePlot(yjsl, features = 'autophagy1',
            label = T, reduction = 'umap')

# Nature journal color standards
{
  nature_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  FeaturePlot(yjsl, features = 'autophagy1', label = TRUE,
              reduction = 'umap', col = nature_colors)
}

VlnPlot(yjsl, features = 'autophagy1')

## Enhance the plot appearance
{
  color <- c('#5470c6', '#91cc75', '#fac858', '#ee6666', '#73c0de', '#3ba272',
             '#fc8452', '#9a60b4', '#ea7ccc')
}

#my_comparisons <- list(c('T cells','B cells'),c('T cells','Mast cells'))
library(ggpubr)
ggviolin(yjsl@meta.data, x = "cellType", y = "autophagy1",
         color = "cellType", add = 'mean_sd', fill = 'cellType',
         add.params = list(color = "black")) +
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend() + labs(x = '')


VlnPlot(yjsl, features = 'autophagy1', col = color, pt.size = 0) #,pt.size = 0


# Divide cells into groups based on the score
yjsl@meta.data$geneSet.Type <- ifelse(
  yjsl@meta.data$autophagy1 > median(yjsl@meta.data$autophagy1),
  "High", "Low"
)
Idents(yjsl) = yjsl@meta.data$geneSet.Type

## Display the groups using two colors
# 1:
{
  DimPlot(yjsl, reduction = 'umap',
          label = F, pt.size = 0.5)
}

# 2:
{
  # Assume the following are the Nature journal color standards
  nature_colors <- c("#0072B2", "#D55E00", "#CC79A7")
  
  DimPlot(yjsl, reduction = 'umap', label = FALSE,
          pt.size = 0.2, group.by = 'geneSet.Type',
          cols = nature_colors)
}


# Identify differentially expressed genes between the high-score and low-score groups
DEGs = FindMarkers(yjsl, ident.1 = 'High', logfc.threshold = 0.3)

write.csv(DEGs, "differentially_expressed_genes_high_vs_low.csv")

## Apply additional filtering if needed
DEGs.filter <- subset(DEGs, p_val_adj < 0.05)
save(yjsl, file = "high_low_scores.rda")

load("high_low_scores.rda")

# Use a violin plot to display a highly variable gene in the high-score and low-score groups
VlnPlot(yjsl, features = 'PRDX1', pt.size = 0) #,pt.size = 0