setwd("")
library(WGCNA)

datExpr0 = read.csv("expression_profile.csv", row.names = 1)

datExpr0 = as.data.frame(t(datExpr0))
datExpr0[1:4,1:4]
dim(datExpr0)


## To reduce the computational workload, select genes with variance in the top 25%
## This step is optional

m.vars = apply(datExpr0, 2, var)
expro.upper = datExpr0[,which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr0 <- data.matrix(expro.upper)
dim(datExpr0)

# Mainly check for missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK
# Continue if TRUE is returned
if (!gsg$allOK){
  # Print the genes or samples containing missing values
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove those with missing values
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


## Sample filtering

sampleTree = hclust(dist(datExpr0), method = "average")

pdf(file = "1_sample_clustering_tree.pdf", width = 10, height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# Remove outliers if present, and set the cutHeight parameter based on the clustering plot

## Separate plot
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  # Replace the values in "h = 110" and "cutHeight = 110" with your desired cutoff
  abline(h = 175, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 175, minSize = 10)
keepSamples = (clust==1)    # Cluster 1 is retained
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

# No removal is needed if there are no outlier samples
datExpr = datExpr0


# Import trait information

# The values must be numeric, such as age, or encoded as 0 and 1, or 1, 2, 3, etc.
## If a sample was removed in the previous step, remember to remove it here as well!!!

datTraits = read.csv("sample_information.csv", row.names = 1)


sampleTree2 = hclust(dist(datExpr), method = "average")
# Trait values for each sample. White indicates low values, red indicates high values, and gray indicates missing values
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot sample clustering and traits together
pdf(file = "2_sample_traits.pdf", width = 10, height = 8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

## Begin the formal WGCNA analysis

## Select the soft threshold
## Set a series of soft thresholds ranging from 1 to 30
powers = c(1:10, seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# This result is the recommended soft threshold
sft$powerEstimat


cex1 = 0.9 # Usually 0.9. It should not be lower than 0.85
pdf(file = "3_soft_threshold_selection.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=cex1, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()


# Proceed with network construction

power = sft$powerEstimate
power # If no soft threshold was recommended above, select one manually based on the plots

## This step takes time

net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned",
                       maxBlockSize = ncol(datExpr) + 1,
                       minModuleSize = 50, # The default minModuleSize is 30. This parameter sets the minimum number of genes per module. Smaller values allow smaller modules to be retained
                       reassignThreshold = 0,
                       mergeCutHeight = 0.6, # The default mergeCutHeight is 0.25. This parameter sets the distance threshold for merging similar modules. Smaller values make merging less likely, retaining more modules
                       deepSplit = 1, # The default deepSplit is 2. This parameter adjusts the sensitivity of module detection. Larger values increase sensitivity and produce more modules
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)


# Display the number of modules identified and the number of genes in each module
table(net$colors)
## If the results are unreasonable, adjust the parameters as appropriate

# For further details, refer to https://zhuanlan.zhihu.com/p/34697561


mergedColors = labels2colors(net$colors)
pdf(file = "4_clustering_trends.pdf", width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Each color represents a gene module, meaning that the genes within it have similar functions
# Gray has no significance

# Save the genes corresponding to each module
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm = data.frame(net$colors)
gm$color = moduleColors
head(gm)

genes = split(rownames(gm), gm$color)
save(genes, file = "genes.Rdata")


# Correlations between modules and traits

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Heatmap
pdf(file = "5_correlation_heatmap.pdf", width = 8, height = 8)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


# The correlation coefficient should preferably exceed 0.8. If none does, 0.6 or 0.7 is also acceptable. Lower values are unsuitable


# Export gene modules to files
text <- unique(moduleColors)
for (i in 1:length(text)) {
  y = t(assign(paste(text[i], "expr", sep = "."), datExpr[moduleColors==text[i]]))
  write.csv(y, paste(text[i], "csv", sep = "."), quote = F)
}


# GS and MM
# GS represents the correlation between each gene in a module and the trait
# MM represents the correlation between each gene and its module, indicating whether it follows the module's trend

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# Set i below to the column number of the trait of greatest interest
# Assign the module most strongly correlated with the trait of interest to module below
traitData = datTraits
i = 2 # Replace with the column number of your desired trait
module = "black" ## Replace with your own sample color
assign(colnames(traitData)[i], traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(instrait), sep="")
names(GSPvalue) = paste("p.GS.", names(instrait), sep="")
pdf(file = paste0("6_module_correlation_scatterplot.pdf"), width = 5, height = 5)
column = match(module, modNames) # Find the column containing the target module
moduleGenes = moduleColors==module # Find the rows corresponding to genes in the module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

## Genes from the core module can also be extracted separately for GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)

hub_gene <- read.csv('black.csv', header = F) ## Replace with your own core module
hub_gene <- hub_gene[,1]
hub_gene <- hub_gene[2:277] # Replace with your own number
gene = unique(hub_gene)

# Convert to ENTREZID
entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
# Use entrezIDs
gene <- entrezIDs$ENTREZID

## GO enrichment analysis
go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="all", readable = T)

write.table(go, file="GO.txt", sep="\t", quote=F, row.names = F) #

## Visualization
## Bar plot
pdf(file="GO-barplot.pdf", width = 10, height = 15)
barplot(go, drop = TRUE, showCategory = 10, label_format=100, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

## Bubble plot
pdf(file="GO-bubble.pdf", width = 10, height = 15)
dotplot(go, showCategory = 10, label_format=100, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

# KEGG analysis
kk <- enrichKEGG(gene = gene, keyType = "kegg", organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "fdr")
write.table(kk, file="KEGG.txt", sep="\t", quote=F, row.names = F)

## Visualization
## Bar plot
pdf(file="KEGG-barplot.pdf", width = 10, height = 13)
barplot(kk, drop = TRUE, showCategory = 15, label_format=100)
dev.off()

## Bubble plot
pdf(file="KEGG-bubble.pdf", width = 10, height = 13)
dotplot(kk, showCategory = 15, label_format=100)
dev.off()


# Genes with high values of both GS and MM can also be identified here
## You can explore this if interested. I have generally not seen these data used in articles

f = data.frame(GS = abs(geneModuleMembership[moduleGenes, column]),
               MM = abs(geneTraitSignificance[moduleGenes, 1]))
rownames(f) = rownames(gm[moduleGenes,])
head(f)


## Display the weighted network using a gene correlation heatmap

# Each row and column represents a gene. Selecting 400 genes is generally sufficient for plotting, otherwise the computer may be unable to handle it

nSelect = 400
set.seed(10)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)

select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# Recalculate the distance-based gene dendrogram. Reclustering is required for a subset of genes
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
library(gplots)
myheatcol = colorpanel(250, 'red', "orange", 'lemonchiffon')
pdf(file = "7_weighted_network.pdf", width = 15, height = 8)
plotDiss = selectTOM^7
diag(plotDiss) = NA # Set the diagonal to NA so that it appears as white points in the plot, making the trends clearer
TOMplot(plotDiss, selectTree, selectColors, col=myheatcol, main = "Network heatmap plot, selected genes")
dev.off()


# Correlations between modules and traits
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(cbind(MEs, instrait))

pdf(file = "8_module_trait_correlations.pdf", width = 12, height = 8)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# The plots above can also be drawn separately

pdf(file = "9_module_clustering_tree.pdf", width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

pdf(file = "10_trait_correlation_heatmap.pdf", width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(4,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()