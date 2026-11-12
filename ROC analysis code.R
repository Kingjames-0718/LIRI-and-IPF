# Load the pROC package for plotting ROC curves
library(pROC)

# Set the working directory


# Read the CSV file, using the first column as row names
rt = read.csv("input_file_single_gene.csv", row.names = 1)

# Select the target variable y, usually the sample class label
y = colnames(rt)[1]
# Change this to your gene of interest
x = "ACTB"

# Generate the ROC curve using the target variable y and the expression values of gene x
ROC = roc(rt[,y], as.vector(rt[,x]))

pdf("ROC_single_gene.pdf", width=5, height=5)

# Plot the ROC curve in red and display the AUC value
plot(ROC, print.auc=TRUE, col="red")
# Close the graphics device to complete the plot
dev.off()