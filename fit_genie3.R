# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(randomForest)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(42, kind = "L'Ecuyer-CMRG")

# Import data
mat <- as.matrix(fread('./dream5/e_coli/net3_expression_data.tsv'))

# Scale (per authors)
mat <- scale(mat)
x <- mat[, seq_len(334)]
y <- mat[, 335:ncol(mat)]

# Loop
rf_loop <- function(gene) {
  f <- randomForest(x, y[, gene], ntree = 1000, mtry = floor(sqrt(334)),
                    importance = TRUE)
  saveRDS(f, paste0('./grn/models/G', 334 + gene, '.rds'))
  out <- f$importance[, 2]
  return(out)
}
imp <- foreach(g = seq_len(ncol(y)), .combine = rbind) %dopar%
  rf_loop(g)
fwrite(imp, './grn/adj_mat.csv')