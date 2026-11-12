# Common gene identifiers in expression matrices
# 1. Gene symbol
# 2. Ensembl Gene ID
# 3. Entrez ID


############## Code starts here ##############

a <- read.table("expression_matrix.txt", sep = "\t", check.names = F, header = T)

library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

a2 <- bitr(a$`Symbol`, fromType = 'SYMBOL',
           toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = "org.Hs.eg.db")

write.csv(a2, "processed_matrix.csv")