# Load the required packages
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)

# Read the data
data <- read.csv("input_data.csv")

# Extract genes from the two columns as sets
comgene <- list(
  'Group1' = unique(data$Group1),
  'Group2' = unique(data$Group2)
)


# Calculate the intersecting genes
com_genes <- intersect(comgene$Group1, comgene$Group2)

# Save as CSV
{
  write.csv(data.frame(Intersection_Genes = com_genes),
            "intersection_genes.csv", row.names = FALSE)
}

# Standard Venn diagram
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # Remove the borders
  fill_color = c("#C2E0F7", "#A4CBA8"),
  set_name_color = c("#C2E0F7", "#A4CBA8"),
  text_color = "black",  # Adjust the text color
  text_size = 5          # Increase the text size
)
print(p1)

# Draw an area-proportional diagram
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.5),  # Adjust the font size and color
  edges = NULL,  # Remove the borders
  fills = c("#C2E0F7", "#A4CBA8"),
  alpha = 0.9,
  quantities = list(cex = 1.5, col = 'black')  # Make the count labels more prominent
)