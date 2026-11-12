# Note: Ensure the folder name matches your local directory
setwd("")

# Load packages required for GSEA analysis
library(clusterProfiler)  # The clusterProfiler package is used for functional enrichment analysis, including GSEA
library(cowplot)          # The cowplot package is used for plot arrangement and combination
library(ggplot2)          # The ggplot2 package is used to draw high-quality visualizations
library(enrichplot)       # The enrichplot package is used to visualize the results of functional enrichment analysis

# Read the input gene expression data file, where the first column serves as row names (i.e., gene names)
rt <- read.csv("Input_file.csv", row.names = 1)

# Load gene sets
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
GSEA <- read.gmt("c5.go.v2024.1.Hs.symbols.gmt") ## Replace with your gene set of interest
# View the first few rows of the gene set data
head(GSEA)


# Sort genes for subsequent GSEA analysis
{
  geneList <- rt$logFC  # Extract the logFC column, which represents the gene expression change value
  names(geneList) <- rownames(rt)  # Set gene names as the names of geneList
  geneList <- sort(geneList, decreasing = T)  # Sort genes from highest to lowest based on expression change values
}

# Perform GSEA enrichment analysis
{
  GSEA_analysis <- GSEA(geneList,             # Input the sorted gene list
                        TERM2GENE = GSEA,     # Gene sets to use
                        pvalueCutoff = 0.05,  # P-value threshold, representing the criteria for significance
                        minGSSize = 20,       # Minimum gene set size in enrichment analysis
                        maxGSSize = 1000,     # Maximum gene set size in enrichment analysis
                        eps = 0,              # Boundary parameter for P-value calculation, usually 0
                        pAdjustMethod = "BH") # P-value correction method, Benjamini-Hochberg method
}
# Convert the GSEA analysis results into a data frame
result <- data.frame(GSEA_analysis)
# View the dimensions of the results
dim(GSEA_analysis@result)  

# Save the analysis results as a CSV file
write.csv(result, "Enrichment_results.csv")


#### The following is the visualization section!!


# Select a single pathway
p1 = gseaplot2(GSEA_analysis, "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_A_NUCLEIC_ACID", 
               color = "red3", # Set curve color
               pvalue_table = T) # Display P-value table in the plot

p1
ggsave(p1, filename = 'Single_gene_pathway.pdf', width = 13, height = 10)

# Draw GSEA enrichment curves for specific pathways —— Multiple pathways
p2 = gseaplot2(GSEA_analysis, c("GOBP_DOUBLE_STRAND_BREAK_REPAIR", 
                                "GOMF_ATP_HYDROLYSIS_ACTIVITY"), 
               color = c("red3", "blue4"), # Set curve colors for different pathways
               pvalue_table = T)

p2

ggsave(p2, filename = 'Multiple_gene_pathways.pdf', width = 13, height = 10)

# Draw a dot plot of enrichment results, displaying the top 10 enriched pathways
p3 = dotplot(GSEA_analysis, showCategory = 10, color = "p.adjust") # Color by adjusted P-value
p3

ggsave(p3, filename = 'Enrichment_pathway_dotplot.pdf', width = 8, height = 8)

# Display activated and suppressed pathways separately, using different colors to distinguish them
{
  p4 = dotplot(GSEA_analysis, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
    theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
          axis.title = element_text(size = 10, color = "black"), 
          axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 0, hjust = 1 ),
          legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
}

p4

ggsave(p4, filename = 'Activation_and_suppression_plot.pdf', width = 10, height = 8)

# Draw a ridge plot to show the distribution of the top 10 enriched pathways
{
  p5 <- ridgeplot(GSEA_analysis,
                  showCategory = 10,
                  fill = "p.adjust",
                  core_enrichment = TRUE,
                  label_format = 30, # Set the character length per line for axis labels; it will wrap automatically if too long.
                  orderBy = "NES",
                  decreasing = F) 
}
p5

ggsave(p5, filename = 'Ridge_plot.pdf', width = 10, height = 8)

# Draw the pathway-gene network plot
{
  p6 = cnetplot(GSEA_analysis, showCategory = 5, # Number of pathways
                foldChange = geneList, 
                colorEdge = "TRUE") # Display different colors for different pathways
}
p6
ggsave(p6, filename = 'Pathway_gene_network.pdf', width = 10, height = 8)

# Calculate the similarity between pathways and draw an enrichment map
net <- pairwise_termsim(GSEA_analysis)

{
  p7 = emapplot(net, showCategory = 10, # showCategory displays the number of items
                color = "pvalue", 
                cex_category = 1, # cex_category adjusts node size
                layout = "kk")
}
p7
ggsave(p7, filename = 'Pathway_correlation.pdf', width = 10, height = 8)

# Rank and plot the top 3 pathways
fun <- lapply(1:3, function(i) {
  anno <- GSEA_analysis[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse = "\n")
  
  gsearank(GSEA_analysis,
           GSEA_analysis$ID[1], 
           title = GSEA_analysis$Description[i]) + 
    xlab(NULL) +
    # ylab(NULL) +
    annotate("text", 5000,
             GSEA_analysis[i, "enrichmentScore"] * .75, 
             label = lab, 
             hjust = 0, vjust = 0)
})
p8 <- plot_grid(plotlist = fun, ncol = 1)

p8

ggsave(p8, filename = 'Pathway_summary_plot.pdf', width = 10, height = 8)