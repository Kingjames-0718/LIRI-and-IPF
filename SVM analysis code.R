## Set the working directory

## Install the corresponding packages
library(tidyverse)
library(glmnet)
source('msvmRFE.R')
library(VennDiagram)
#BiocManager::install("sigFeature")
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)


## Read the file
input <- read.csv("processed_data.csv", row.names = 1)

## Begin the formal analysis
### Use five-fold cross-validation
svmRFE(input, k = 5, halve.above = 100)

nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)

# # View the main variables
top.features = WriteFeatures(results, input, save=F)
head(top.features)

# Save the features identified by SVM-RFE to a file
write.csv(top.features, "feature_svm.csv")

## The more variables, the longer the execution time
featsweep = lapply(1:20, FeatSweep.wrap, results, input)

save(featsweep, file = "result.RData")
load("result.RData")

# Plot the results
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

# Plot the error rate
PlotErrors(errors, no.info=no.info)

# Plot the accuracy
Plotaccuracy(1-errors, no.info=no.info)


# The position marked by the red circle in the plot indicates the lowest error rate
which.min(errors)

## Save the core gene set
top <- top.features[1:which.min(errors), "FeatureName"]

write.csv(top, "selected_core_genes.csv")