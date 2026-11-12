# Note: Ensure the folder name matches your local directory
setwd("")

# Import required R packages
library(dplyr)          # Data manipulation package
library(CIBERSORT)      # For immune cell infiltration analysis
library(ggplot2)        # Draw data charts
library(pheatmap)       # Draw heatmaps
library(ggpubr)         # Provides ggplot2 enhancements, convenient for plotting and statistics
library(reshape2)       
library(tidyverse)      
library(limma)          # For differential expression analysis
library(ggthemes)       
library(RColorBrewer)   # Provides color palettes

# Read expression matrix and group information
expr <- read.csv("Expression_matrix.csv", row.names = 1)  # Read gene expression matrix, first column is gene ID
group <- read.csv("Group.csv", row.names = 1)             # Read group information (e.g., high and low groups)

## CIBERSORT has strict requirements for data format:
# (1) Cannot have negative or missing values; (2) Standardize data as much as possible
boxplot(expr, outline = FALSE, notch = F, las = 2)  # Draw data boxplot to check data distribution
expr = normalizeBetweenArrays(expr)  # Data standardization to eliminate batch effects
boxplot(expr, outline = FALSE, notch = F, las = 2)  # Boxplot of data after standardization

# Immune cell types (downloaded from the official website)
LM22 <- read.table("LM22.txt", header = TRUE, row.names = 1, sep = "\t") %>% as.matrix()

# Start CIBERSORT analysis, which takes a relatively long time
# The calculated results contain the abundance of 22 immune cells, and the last three columns are other statistics, which can be ignored.
# perm: represents the number of permutations, the larger the number, the longer the running time, generally set to 1000
# QN: set to TRUE for microarray data, FALSE for sequencing data
result <- cibersort(sig_matrix = LM22, mixture_file = expr, perm = 1000, QN = FALSE) %>% as.data.frame()
result1 <- result 
# Save the CIBERSORT analysis results to a CSV file
write.csv(result, "CIBERSORT_results.csv")

# Keep only the abundance data of 22 immune cell types, remove the last three columns of statistical information
result <- result[, 1:22]

# Draw a scatter plot of immune cell type correlations
# This can also be placed at the end to only plot immune cells with statistical significance
library(corrplot)  
cor <- cor(result)  # Calculate the correlation of immune cell abundance
corrplot(cor, method = 'number')  # Display correlation values
# Another layout
corrplot(cor, order = "AOE", type = "upper", tl.pos = "d")



# Merge immune cell abundance data with group information
data1 <- cbind(result, group)  # Merge immune cell abundance and group information
data1 <- data1 %>% rownames_to_column("sample")  # Convert row names (sample IDs) to a column
data <- gather(data1, key = CIBERSORT, value = Proportion, -c(group, sample))  # Convert data format for easier plotting

## Draw boxplots for group differences

{
  ggboxplot(data, x = "CIBERSORT", y = "Proportion",  # Draw boxplot of immune cell proportions
            fill = "group", palette = c("#028A8B", "#FB6346")) +  # Set group colors
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +  # Rotate x-axis labels
    stat_compare_means(label = "p.signif", method = "t.test",  # Perform t-test and annotate p-values
                       ref.group = ".all.", hide.ns = TRUE)  # Compare differences among all groups
}

# Draw heatmap of immune cell abundance
{
  k <- apply(result, 2, function(x) { sum(x == 0) < nrow(result1) / 2 })  # Filter out non-zero immune cell types
  re2 <- as.data.frame(t(result[, k]))  # Transpose data for heatmap plotting
  an = data.frame(group = group, row.names = colnames(expr))  # Create column annotation data
}

pheatmap(re2, scale = "row", show_colnames = FALSE, cluster_cols = FALSE,  # Draw heatmap
         annotation_col = an, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))  # Set heatmap colors

# Extract differences for a single immune cell (e.g., CD8 T cells)
dat2 = data[data$CIBERSORT == "T.cells.CD8",]  # Extract proportion data for CD8 T cells

{
  ggplot(dat2, aes(x = group, y = Proportion)) +  # Draw boxplot for CD8 T cells
    labs(y = "Cell composition", x = NULL, title = "") +
    geom_boxplot(aes(fill = group), position = position_dodge(0.5), width = 0.5, size = 0.4,
                 outlier.alpha = 1, outlier.size = 0.5) +  # Set boxplot style
    theme_bw() + 
    scale_fill_manual(values = c("#028A8B", "#FB6346")) +  # Set group colors
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test",  # Use Wilcoxon test
                       size = 6, hide.ns = TRUE)
}

# Draw stacked proportion plot to show the proportion of immune cells in each sample
{
  identical(rownames(result), group$Samples)  # Confirm that immune cell data and sample group information match
  data <- cbind(rownames(result), result)  # Merge sample names and immune cell abundance data
  colnames(data)[1] <- "Samples"  # Rename the first column after merging to "Samples"
  data <- melt(data, id.vars = c("Samples"))  # Convert data from wide format to long format
  colnames(data) <- c('Samples', 'celltype', 'proportion')  # Set column names
}

# Draw stacked proportion plot
{
  colour = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"))  # Set colors
  ggplot(data, aes(Samples, proportion, fill = celltype)) +  # Draw stacked bar chart
    geom_bar(stat = "identity", position = "fill") +  # Draw stacked proportion plot
    scale_fill_manual(values = colour) +  # Set colors
    ggtitle("") + theme_gray() + theme(axis.ticks.length = unit(3, 'mm'), axis.title.x = element_text(size = 11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    guides(fill = guide_legend(title = ""))
}

# Another representation of the stacked plot
{
  ggplot(data, aes(x = Samples, y = proportion, fill = celltype)) +  # Draw stacked bar chart
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = colour) +  # Set colors
    labs(x = "", y = "", title = "") +  # Set title and axis labels
    scale_y_continuous(expand = c(0, 0)) +  # Remove white space on y-axis
    guides(fill = guide_legend(ncol = 1)) +  # Set legend arrangement
    theme_bw() + theme(legend.key = element_blank(), legend.title = element_blank(),
                       panel.grid = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5))
}