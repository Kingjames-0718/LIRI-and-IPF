library(magrittr)
library(reshape2)
library(tidyr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(cowplot)

# Import the data
load("yjsl_Type.rda")

# Check the cell types
table(yjsl@meta.data$cellType, yjsl@meta.data$orig.ident)

# Create a table of samples, cell types, and cell counts
df <- table(yjsl@meta.data$cellType, yjsl@meta.data$orig.ident) %>% melt()

# Customize the column names
colnames(df) <- c("Cluster", "Sample", "Number")

# Define the cell type names
cluster = c("T cells", "Macrophage cells", "Epithelial cells",
            "B cells", "Fibroblast cells", "Mast cells",
            "Endothelial cells")

df$Cluster <- factor(df$Cluster, levels = cluster)

# Define colors
{
  color <- c("#d95f02", "#66a61e", "#1b9e77", "#e7298a", "#386cb0", "red3", "yellow4")
  df$Sample <- factor(df$Sample) # Convert the Sample column to a factor
}

## Proportions of samples within each cell type
{
  p1 <- ggplot(data = df, aes(x = Number, y = Cluster, fill = Sample)) +
    geom_bar(stat = "identity", width = 0.8, position = "fill") +
    scale_fill_manual(values = color) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "", y = "Ratio") +
    theme(axis.text.y = element_text(size = 12, colour = "black")) +
    theme(axis.text.x = element_text(size = 12, colour = "black"))
}
p1

# Proportions of cell types within each sample
{
  p2 <- ggplot(data = df, aes(x = Number, y = Sample, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.8, position = "fill") +
    scale_fill_manual(values = color) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "", y = "Ratio") +
    theme(axis.text.y = element_text(size = 12, colour = "black")) +
    theme(axis.text.x = element_text(size = 12, colour = "black"))
}
p2

## Pie chart
{
  names <- table(yjsl$cellType) %>% names()
  ratio <- table(yjsl$cellType) %>% as.numeric()
  pielabel <- paste0(names, " (", round(ratio / sum(ratio) * 100, 2), "%)")
}

pdf(file = "pie_chart.pdf", width = 10, height = 10)

pie(ratio, labels = pielabel,
    radius = 1.0, clockwise = T,
    main = "cellType", col = color)

dev.off()