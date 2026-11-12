setwd("")


library(AUCell)
library(ggplot2)
library(Seurat)
library(clusterProfiler)


load("yjsl_Type.rda")

DimPlot(yjsl, group.by = 'cellType', label = TRUE)

# Randomly sample cells
sc.id = sample(colnames(yjsl), 3000)
yjsl = yjsl[, sc.id]

# Rank genes
cells_rankings <- AUCell_buildRankings(
  exprMat = yjsl@assays$RNA@data,
  featureType = "genes",
  splitByBlocks = TRUE,
  nCores = 2,
  plotStats = TRUE
)

# Read the gene sets
c2 <- read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
geneSets <- lapply(unique(c2$term), function(x){print(x); c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)

# Calculate AUC values
cells_AUC <- AUCell_calcAUC(
  geneSets,
  cells_rankings, nCores = 1,
  aucMaxRank = nrow(cells_rankings) * 0.1
)

pdf("AUC histogram.pdf")
cells_assignment <- AUCell_exploreThresholds(
  cells_AUC,
  plotHist = TRUE, assign = TRUE
)

dev.off()

# Select a pathway of interest for visualization
library(msigdbr)

## Filter by keyword
#grep("RENIN", rownames(cells_AUC@assays@data$AUC), value = T)

geneSet <- "KEGG_RENIN_ANGIOTENSIN_SYSTEM"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
yjsl$AUC <- aucs


## Plot the results
# 1
{
  plot.df <- data.frame(
    yjsl@meta.data,
    yjsl@reductions$umap@cell.embeddings
  )
  
  ggplot() +
    geom_point(
      data = plot.df,
      aes(x = UMAP_1, y = UMAP_2, colour = AUC),
      size = 0.5
    ) +
    viridis::scale_color_viridis(option = "A") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(colour = 'black', size = 14),
      axis.title = element_text(size = 14),
      axis.line = element_line(colour = 'black'),
      panel.border = element_rect(
        size = 0.5, linetype = "solid", colour = "black"
      ),
      aspect.ratio = 1
    )
}

DimPlot(yjsl, group.by = 'cellType', label = TRUE)

# 2
{
  class_avg <- plot.df %>%
    group_by(cellType) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  ggplot(plot.df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = AUC)) +
    viridis::scale_color_viridis(option = "A") +
    ggrepel::geom_label_repel(
      aes(label = cellType),
      data = class_avg,
      size = 6,
      label.size = 0,
      segment.color = NA
    ) +
    theme(legend.position = "none") +
    theme_bw()
}


highcelltype <- c("Mast cells")

# Divide into high and low groups
{
  yjsl@meta.data$geneSet.Type <- ifelse(
    yjsl@meta.data$cellType %in% highcelltype, "High", "Low"
  )
}

DimPlot(
  yjsl, reduction = "umap", group.by = "geneSet.Type",
  label = TRUE, repel = TRUE
)

{
  Idents(yjsl) <- "geneSet.Type"
  DEGs <- FindAllMarkers(
    yjsl,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

write.csv(DEGs, "differentially_expressed_genes_high_vs_low.csv")