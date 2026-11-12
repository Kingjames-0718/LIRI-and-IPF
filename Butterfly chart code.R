# install.packages("devtools")
# devtools::install_github("Hy4m/linkET", force = TRUE)
# Load the required R packages
library("linkET")
library(ggplot2)
library(vegan)

# Read the expression data
# file.choose() opens a file selection dialog for manually selecting a file
varechem <- read.delim(file.choose(), head = T, row.names = NULL)

# Read the gene expression/risk file data
varespec <- read.delim(file.choose(), head = T, row.names = NULL)

# Calculate the Pearson correlation matrix between the datasets
{
  correlate(varechem, varespec) %>%
    # Use qcorrplot to draw a correlation heatmap
    qcorrplot() +
    # Add square cells to improve the appearance of the heatmap
    geom_square() +
    # Set the color gradient using the RColorBrewer "RdBu" palette with 11 colors
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
}

# Draw only the lower triangle of the correlation heatmap to avoid duplicate display
qcorrplot(varechem, type = "lower") +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))


library(dplyr)  # Load the data manipulation package
# Perform the Mantel test to analyze the correlation between expression and risk factors
mantel <- mantel_test(varespec, varechem, spec_select = list(Riskscore = 1)) %>%
  # Categorize the Mantel test results into intervals
  mutate(
    # Divide the correlation coefficient r into three intervals of <0.2, 0.2-0.4, and >=0.4
    rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    # Divide p-values into three intervals of <0.01, 0.01-0.05, and >=0.05
    pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )


# Draw the correlation heatmap with the lower triangle only and the diagonal hidden
{
  qcorrplot(correlate(varechem), type = "lower", diag = FALSE) +
    geom_square() +
    # Add connections representing the Mantel test results, with color indicating p-values and thickness indicating r-values
    geom_couple(
      aes(colour = pd, size = rd),  # Map color to p-value intervals and size to r-value intervals
      data = mantel,                # Use the Mantel test results
      curvature = nice_curvature()  # Set the curvature of the connections
    ) +
    # Set the heatmap color gradient
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
    # Set the mapped values for connection thickness
    scale_size_manual(values = c(0.5, 1, 2)) +
    # Set the mapped values for connection colors
    scale_colour_manual(values = color_pal(3)) +
    # Customize the legend styles and titles
    guides(
      size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
      colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
      fill = guide_colorbar(title = "Pearson's r", order = 3)
    )
  }