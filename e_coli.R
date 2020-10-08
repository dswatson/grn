# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(randomForest)
library(kernlab)
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
out <- foreach(g = seq_len(ncol(y)), .combine = rbind) %dopar%
  rf_loop(g)
fwrite(out, './grn/adj_mat.csv')

### BASELINE SIMULATION ###

# Simulate baseline x's
Sigma <- cov(x)
sim_x_fn <- function(n) {
  mu <- matrix(rnorm(n * 334), nrow = n)
  sim_x <- mu %*% chol(Sigma)
  colnames(sim_x) <- colnames(x)
  return(sim_x)
}
sim_x <- sim_x_fn(n = 12000)

# Simulate baseline y's
y_loop <- function(sim_x, y_gene) {
  f <- readRDS(paste0('./grn/models/G', 334 + y_gene, '.rds'))
  y_hat <- predict(f, sim_x)
  return(y_hat)
}
sim_y <- foreach(g = seq_len(ncol(y)), .combine = cbind) %dopar%
  y_loop(sim_x, g)
colnames(sim_y) <- colnames(y)

# Output baseline data
sim_pre <- cbind(sim_x, sim_y)

### INTERVENTION SIMULATION ###

# Simulate the effect of random knockouts
ko_loop <- function(n, x_gene) {
  sim_x_prime <- sim_x_fn(n)
  sim_x_prime[, x_gene] <- min(x[, x_gene]) - 1
  sim_y <- foreach(g = seq_len(ncol(y)), .combine = cbind) %dopar%
    y_loop(sim_x_prime, g)
  colnames(sim_y) <- colnames(y)
  sim_post <- cbind(sim_x_prime, sim_y)
  fwrite(sim_post, paste0('./grn/simulations/ko_G', x_gene, '.csv'))
}

# Compute phis
adj_mat <- ifelse(out >= 10, 1, 0)
outdegree <- colSums(adj_mat)
keep <- outdegree > 100
phi_fn <- function(sim_x, sim_y, tf) {
  x <- sim_x[, tf]
  y <- sim_y[, adj_mat[, tf] == 1]
  dat <- cbind(x, y)
  kf <- rbfdot(sigma = 1e-4) # EH??
  k_mat <- kernelMatrix(kernel = kf, x = dat)
  pca <- kpca(k_mat) 
  out <- rotated(pca)[, 1]
  return(out)
}







q <- sapply(seq_len(nrow(a)), function(i) {
  quantile(a[i, ], 0.95)
})

b <- ifelse(a >= 10, 1, 0)





# Simulate baseline orphans
indegree <- rowSums(a)
emp_orphans <- y[, indegree == 0]
Sigma <- cov(emp_orphans)
mu <- matrix(rnorm(n * ncol(emp_orphans)), nrow = n)
sim_orphans <- mu %*% chol(Sigma)
colnames(sim_orphans) <- colnames(emp_orphans)







