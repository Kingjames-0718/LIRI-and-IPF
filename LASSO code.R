library(tidyverse)
library(broom)
library(glmnet)


# Read the expression matrix and classification file, using the first column as row names
x <- read.csv("expression_matrix.csv", row.names = 1)
y <- read.csv("classification_file.csv", row.names = 1)

# Convert the imported data to matrices for subsequent model calculations
x = as.matrix(x)
y = as.matrix(y)

# Set the random seed to ensure reproducibility
set.seed(66666)
#?cv.glmnet

# Use cross-validation to select the optimal LASSO model parameter
# type.measure = "mse" specifies mean squared error as the model evaluation metric
# nfolds = 5 specifies five-fold cross-validation
# alpha = 1 specifies LASSO regression, while alpha = 0 specifies Ridge regression
cvfit = cv.glmnet(x, y, type.measure = "mse", nfolds = 5, alpha = 1)

# Plot the cross-validation results to show model performance at different lambda values
plot(cvfit)


# View the minimum lambda value calculated during cross-validation, corresponding to the optimal model
cvfit$lambda.min

c(cvfit$lambda.min, cvfit$lambda.1se)

# Refit the LASSO model using a wider range of lambda values
# family = "binomial" indicates that the dependent variable is a binary classification variable
lasso <- glmnet(x, y, family = "binomial", alpha = 1, nlambda = 100) # family = "binomial" is suitable for binary discrete dependent variables

# Plot the LASSO regression paths to show how coefficients change at different lambda values
plot(lasso, xvar = "lambda", label = F) # The lambda parameter


# Extract the model coefficients at the minimum lambda value
coef <- coef(lasso, s = cvfit$lambda.min)
# Identify the indices of nonzero coefficients to select genes that contribute to the model
index <- which(coef != 0)
actCoef <- coef[index]
# Extract the gene names and coefficient values corresponding to nonzero coefficients
lassoGene = row.names(coef)[index]
geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
# Save the selected core genes and their coefficients to a CSV file
write.csv(geneCoef, "core_gene_coefficients.csv")