# Load the required R packages
library(GEOquery)
library(stringr)
library(limma)
library(tidyverse)   #

# The getGEO function downloads the dataset with the specified ID from the GEO database. "GSE17351" is the dataset ID
gset <- getGEO("GSE17351", destdir=".", AnnotGPL = F, getGPL = F)

# Extract the first expression matrix from the downloaded dataset
a <- gset[[1]]
# Obtain the expression matrix
dat <- exprs(a)

# Check the dimensions of the data matrix, representing the numbers of genes and samples
dim(dat)
# Display the first four rows and four columns of the expression matrix to help check the data format
dat[1:4, 1:4]


# Manually download the platform data
ids = read.table('annotation.txt', sep = "\t", header = T)

# Check the data format
ids[1:5,]

# Extract the first gene before "///" in the Gene.Symbol column
{
  ids$Gene.Symbol <- sapply(strsplit(as.character(ids$Gene.Symbol), " /// "), `[`, 1)
}
# View the first five rows
head(ids, 5)

{
  # Rename the columns of the annotation data frame to "probe_id" (probe ID) and "symbol" (gene symbol)
  colnames(ids) <- c("probe_id", "symbol")
  
  # View the first few rows
  head(ids)
}

# Match the probe IDs to the row names of the expression matrix
ids <- ids[match(rownames(dat), ids$probe_id), ]
# Match the annotation data frame to the expression matrix by probe ID. The match function returns the matching positions
# This ensures that the row names of the expression matrix (probe IDs) correspond to the correct gene symbols

# Remove duplicate gene symbols. If multiple probes map to the same gene, retain one gene symbol
ids <- ids[!duplicated(ids$symbol), ]


# Reorder the rows of the expression matrix to match the order of probe IDs in ids
dat <- dat[ids$probe_id, ]
rownames(dat) <- ids$symbol  # Replace the original probe IDs with gene symbols as row names

dat <- as.data.frame(dat)  # Convert dat to a data frame to ensure consistency in data operations

# Remove genes with column names of NA
dat = na.omit(dat)
# Use na.omit to remove rows containing NA (missing values) from the expression matrix, thereby removing genes with missing values

# Remove low-expression genes
dat = dat[rowMeans(dat) > 1,]
# Remove genes with an average expression value below 1 across all samples, thereby removing low-expression genes

# If expression values exceed 15, run the following code
dat = log2(dat + 1)
# Apply a log2 transformation to the expression matrix to ensure an appropriate numerical scale. Adding 1 avoids taking the logarithm of zero

# Data correction ------------------------------
boxplot(dat, outline=FALSE, notch=F, las=2)


dat = normalizeBetweenArrays(dat)
# Use the normalizeBetweenArrays function from the limma package to normalize the data and eliminate technical variation between samples

boxplot(dat, outline=FALSE, notch=F, las=2)
# Draw another box plot after normalization to check whether the data distributions are more consistent

# Save as a text file
write.table(dat, file = "GSE17351.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)