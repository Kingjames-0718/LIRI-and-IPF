library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(harmony)
library(plyr)
library(randomcoloR)

## A package specifically for enhancing plot appearance
#install.packages("/Users/zhangzhenhu/Desktop/single_cell_annotation_manual/scRNAtoolVis_0.1.0.tar.gz", repos=NULL, type="source")
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

yjsl = readRDS("af2.rds")

################ Recommended method of manual annotation ################
###### First approach using differentially expressed genes ################
pdf(file = "figure_name.pdf", width = 6, height = 6)

DimPlot(yjsl, reduction = "umap", label = T)

dev.off()

# Identify all differentially expressed genes
Idents(yjsl) <- yjsl$seurat_clusters
allmarkers <- FindAllMarkers(yjsl, only.pos = TRUE,
                             min.pct = 0.25, # Only genes expressed in at least 25% of cells are considered for differential expression
                             logfc.threshold = 0.25) # Select only genes with a log2 fold change greater than 0.25

save(allmarkers, file = "allmarkers.rda")
load("allmarkers.rda")

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

yjsl <- FindClusters(yjsl, resolution = 0.1) ## Remember to use the same resolution value as above

genes <- list(
  "T cells" = c("CD3D", "CD3E", "TRAC"),
  "Monocytes/macrophages" = c("CD163", "CD14"),
  "B cells" = c("CD79A", "MS4A1"),
  "Fibroblasts" = c("FGF7", "MME"),
  "Endothelial cells" = c("PECAM1", "VWF"),
  "mast cells" = c("TPSAB1", "TPSB2"),
  "Epithelial_cells" = c("CDH1", "CLDN4")
)

do_DotPlot(sample = yjsl, features = genes, dot.scale = 10,
           colors.use = c("yellow", "red"), legend.length = 50,
           legend.framewidth = 2, font.size = 10)

# Bubble plot with another color scheme
DotPlot(yjsl, features = genes, cols = "RdYlBu") +
  RotatedAxis()

# Manual annotation
table(yjsl@active.ident)
ann.ids <- c(
  "T cells",             # cluster0
  "Macrophage cells",    # cluster1
  "Epithelial cells",    # Continue in order below
  "B cells",
  "Macrophage cells",
  "Fibroblast cells",
  "Epithelial cells",
  "T cells",
  "Mast cells",
  "Endothelial cells",
  "Fibroblast cells"
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

library(ggsci)
cors <- pal_igv()(12) # Define colors

DimPlot(yjsl, reduction = "tsne", label = T, group.by = "seurat_clusters", cols = cors)
DimPlot(yjsl, reduction = "umap", label = T, group.by = "cellType", cols = cors)

## Draw a Nature-quality UMAP plot
clusterCornerAxes(
  object = yjsl, reduction = 'umap',
  clusterCol = "cellType", # Classification criterion
  #noSplit = F, groupFacet = 'cellType', # Use when grouping
  addCircle = F, cicAlpha = 0.03, nbin = 100, # Draw circles
  cellLabel = T, cellLabelSize = 3.5 # Labels
) +
  scale_color_igv() + scale_fill_igv()

#### Alternatively, use the following approach ####
## Define a color vector containing the color codes to be used for cell type visualization
my.color <- c(
  "#AF4034", "#55967e", '#006a8e', "#6a60a9",
  "#D55E00", "#E39E3E", "#0072B2",
  "#CC79A7", '#f9a11b', "#4F86C6", "#fdc23e",
  "#e3632d"
)

p1 <- DimPlot(
  object = yjsl, cols = my.color, pt.size = 0.5,
  label = TRUE, reduction = "umap", label.size = 4, repel = TRUE
) +
  theme(text = element_text(family = "sans"))
p1

# Explicitly save the plot with an appropriate resolution and figure size
ggsave("DimPlot_umap.pdf", plot = p1, width = 10, height = 8, dpi = 300)


############# Identify differentially expressed genes again #######################

markers_Type <- FindAllMarkers(
  yjsl, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25
)

save(markers_Type, file = "markers_Type.rda")

load("markers_Type.rda")

write.csv(markers_Type, file = "markers_after_clustering.csv")

# Top 5 highly expressed genes
top5markers_Type <- markers_Type %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Heatmap
DoHeatmap(yjsl, features = top5markers_Type$gene) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(
    low = '#0099CC', mid = 'white', high = '#CC0033',
    name = 'Z-score'
  )

####### Define colors ##########
colaa = distinctColorPalette(100)

# Standard violin plot
VlnPlot(yjsl, features = "CD3E", group.by = "cellType") + NoLegend()

## UMAP plot for a single gene
FeaturePlot(yjsl, features = "CD3E")

# Distribution plots for multiple genes
symbol = c(
  "CD3D", # Enter the genes you want to display
  "CD163",
  "CD79A",
  "FGF7",
  "PECAM1",
  "TPSAB1",
  "CDH1"
)

VlnPlot(
  yjsl, features = symbol, group.by = "cellType",
  stack = TRUE, cols = colaa
) + NoLegend()

# Split violin plots by group
p <- VlnPlot(
  yjsl, features = symbol, stack = T, pt.size = 0,
  flip = T, split.by = 'cellType',
  group.by = "cellType",
  cols = c("#78C2C4", "#C73E3A"),
  split.plot = T
) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = 'black', size = 10, angle = 90),
    legend.position = 'none'
  )
p

# Save the file with assigned cell types for subsequent cell-cell communication and pseudotime trajectory analyses
save(yjsl, file = "yjsl_Type.rda")