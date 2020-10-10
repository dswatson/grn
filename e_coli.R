# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(randomForest)
library(kernlab)
library(Rfast)
library(tidyverse)
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

### SIMULATION FUNCTIONS ###

# Simulate x function
Sigma <- cov(x)
sim_x_fn <- function(n) {
  mu <- matrix(rnorm(n * 334), nrow = n)
  sim_x <- mu %*% chol(Sigma)
  colnames(sim_x) <- colnames(x)
  return(sim_x)
}

# Simulate y function
sim_y_fn <- function(sim_x, y_gene) {
  f <- readRDS(paste0('./grn/models/G', 334 + y_gene, '.rds'))
  y_hat <- predict(f, sim_x)
  return(y_hat)
}

# Wrapper function
sim_dat <- function(n, ko) {
  sim_x <- sim_x_fn(n)
  if (ko > 0) {
    sim_x[, ko] <- min(x[, ko]) - 1
  }
  sim_y <- foreach(g = seq_len(ncol(y)), .combine = cbind) %dopar%
    sim_y_fn(sim_x, g)
  colnames(sim_y) <- colnames(y)
  exprs <- cbind(sim_x, sim_y)
  out <- as.data.table(exprs)[, w := paste0('w', ko)]
}
sim_out <- foreach(g = 0:10, .combine = rbind) %do%
  sim_dat(n = 2000, g)
saveRDS(sim_out, './grn/simulations/sim_exprs.rds')

# Compute phis
adj_mat <- ifelse(imp >= 10, 1, 0)
outdegree <- colSums(adj_mat)
keep <- which(outdegree > 100)

# Kernel PCA
phi_fn <- function(dat, tf) {
  dat <- as.data.frame(dat)
  x <- dat[1:2000, tf]
  y <- dat[1:2000, 334 + which(adj_mat[, tf] == 1)]
  trn <- as.matrix(cbind(x, y))
  d <- Dist(trn) %>% keep(lower.tri(.))
  s <- 1 / median(d)
  pca <- kpca(trn, kernel = 'rbfdot', kpar = list(sigma = s), features = 1)
  w0 <- rotated(pca)
  x_tst <- dat[2001:nrow(dat), tf]
  y_tst <- dat[2001:nrow(dat), 334 + which(adj_mat[, tf] == 1)]
  tst <- as.matrix(cbind(x_tst, y_tst))
  w_rest <- predict(pca, tst)
  out <- c(w0, w_rest)
  return(out)
}
phi_mat <- foreach(g = keep, .combine = cbind) %dopar%
  phi_fn(sim_out, g)
colnames(phi_mat) <- paste0('phi', seq_len(length(keep)))

# Export
res <- cbind(sim_out, phi_mat)
fwrite(res, './grn/simulations/sim_dat.csv')






phi_fn <- function(dat, dat2, tf) {
  dat <- as.data.frame(dat)
  x <- dat[1:2000, tf]
  y <- dat[1:2000, 334 + which(adj_mat[, tf] == 1)]
  trn <- as.matrix(cbind(x, y))
  d <- Dist(trn) %>% keep(lower.tri(.))
  s <- 1 / median(d)
  pca <- kpca(trn, kernel = 'rbfdot', kpar = list(sigma = s), features = 1)
  x_tst <- dat2[, tf]
  y_tst <- dat2[, 334 + which(adj_mat[, tf] == 1)]
  tst <- as.matrix(cbind(x_tst, y_tst))
  out <- as.numeric(predict(pca, tst))
  return(out)
}
phi_mat <- foreach(g = keep, .combine = cbind) %dopar%
  phi_fn(sim_out, baseline, g)
colnames(phi_mat) <- paste0('phi', seq_len(length(keep)))

