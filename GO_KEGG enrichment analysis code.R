library(clusterProfiler) # GO enrichment analysis, KEGG pathway enrichment analysis
library(org.Hs.eg.db)    # Gene annotation database
library(enrichplot)
library(ggplot2)
library(GOplot)

# Set working directory
setwd("")

# Read data
genes_df <- read.csv("genes.csv")  

# Convert gene symbols to ENTREZ IDs
entrezIDs <- bitr(genes_df$gene, fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db) # The 'OrgDb' parameter specifies the human genome annotation database to use

# Use entrezIDs 
gene <- entrezIDs$ENTREZID

## GO enrichment analysis
go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", readable = T)
write.table(go, file = "GO_results.txt", sep = "\t", quote = F, row.names = F) 

## Visualization
## Bar plot
pdf(file = "GO_Barplot.pdf", width = 8, height = 7)
## 'showCategory' changes the number of items you want to display
barplot(go, drop = TRUE, showCategory = 6, label_format = 100, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = 'free')
dev.off()

## Dot plot
pdf(file = "GO_Dotplot.pdf", width = 8, height = 7)
dotplot(go, showCategory = 6, label_format = 100, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = 'free')
dev.off()


# KEGG analysis
kk <- enrichKEGG(gene = gene, keyType = "kegg", organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "fdr")   
write.table(kk, file = "KEGG_results.txt", sep = "\t", quote = F, row.names = F)                       

## Visualization
## Bar plot
pdf(file = "KEGG_Barplot.pdf", width = 8, height = 6)
barplot(kk, drop = TRUE, showCategory = 8, label_format = 100)
dev.off()

## Dot plot
pdf(file = "KEGG_Dotplot.pdf", width = 8, height = 6)
dotplot(kk, showCategory = 8, label_format = 100)
dev.off()