## limma Differential Expression Analysis and Volcano Plot Code ##

# Note: Ensure the folder name matches your local directory
setwd("")

# Load required packages
library(tidyverse)
library(limma)

exp <- read.csv("Expression_matrix_new.csv", row.names = 1) # Read the expression matrix file

## Data cleaning ###
## Not mandatory to run, run only if necessary ###

# If expression values exceed 15, you may need to run the following code for log transformation ------------------------------
exp = log2(exp + 1)  

# Standardize if means are inconsistent ------------------------------
boxplot(exp, outline = FALSE, notch = F , las = 2)
exp = normalizeBetweenArrays(exp)
boxplot(exp, outline = FALSE, notch = F , las = 2)

## Filter out low expression genes ------------------------------

exp = exp[rowMeans(exp) > 1, ] # Filter out low-expressed genes based on your needs; other thresholds can also be applied

# Read classification file ------------------------------

group_info <- read.csv("Classification_new.csv", row.names = 1) # Read the classification file corresponding to the expression matrix headers

## Start differential expression analysis -----------------------------

group_list <- factor(group_info$group, levels = c("Normal", "Tumor"))
design <- model.matrix(~group_list)
con <- lmFit(exp, design)
con2 <- eBayes(con)
DEG1 <- topTable(con2, coef = 2, number = Inf)
DEG2 = na.omit(DEG1) 

## Export differential analysis results -----------------------------

write.csv(DEG2, "Differential_analysis_results_new.csv")


# Filter differentially expressed genes ----------------------------

logFC_cut = 1
p_cut = 0.05

type1 = (DEG2$P.Value < p_cut) & (DEG2$logFC < -logFC_cut)
type2 = (DEG2$P.Value < p_cut) & (DEG2$logFC > logFC_cut)

DEG2$type = ifelse(type1, "Down", ifelse(type2, "Up", "NOT"))

head(DEG2)

write.csv(DEG2, "Differential_analysis_annotated_new.csv")


## Draw volcano plot ------------------- with marker genes ---------

library(ggplot2)

# Add a new column to store label information, list the labels you want to display:
#DEG2$label <- NA
#DEG2$label[rownames(DEG2) %in% c("HIST1H1E")] <- "HIST1H1E" # If HIST1H1E exists

ggplot(DEG2, aes(logFC, -log10(P.Value))) +
  # Horizontal reference line:
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  # Vertical reference line:
  geom_vline(xintercept = c(-0, 0), linetype = "dashed", color = "#999999") +
  # Scatter plot:
  geom_point(aes(size = -log10(P.Value), color = -log10(P.Value))) +
  # Specify color gradient mode:
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  # Specify scatter size gradient mode:
  scale_size_continuous(range = c(1, 3)) +
  # Add labels
  #geom_text(aes(label = label), vjust = 1.5, hjust = 0.5) +
  # Theme adjustments:
  theme_bw() +
  theme(panel.grid = element_blank())