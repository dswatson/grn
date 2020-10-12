# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(ranger)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import, prep data
z <- readRDS('./grn/simulations/baseline.rds')
phis <- readRDS('./grn/simulations/phis.rds')

# Adjacency matrix
imp <- fread('./grn/adj_mat.csv')
adj_mat <- ifelse(imp >= 10, 1, 0)
outdegree <- colSums(adj_mat)
keep <- which(outdegree >= 100)

# Random forest loop
phi_loop <- function(phi_idx) {
  trn_idx <- seq_len(8000)
  tf <- keep[phi_idx]
  x <- z[, c(tf, 334 + which(adj_mat[, tf] == 1))] %>%
    as.data.frame(.) %>%
    mutate(w = paste0('w', rep(1:10, each = 1000)))
  p <- ncol(x)
  f <- ranger(x = x[trn_idx, ], y = phis[trn_idx, phi_idx], 
              mtry = p/3, num.trees = 500, num.threads = 8)
  saveRDS(f, paste0('./grn/phi_models/phi', phi_idx, '.rds'))
  y_tst <- phis[-trn_idx, phi_idx]
  y_hat <- predict(f, x[-trn_idx, ])$predictions
  out <- data.table(
    'phi' = paste0('phi', tf),
     'r2' = cor(y_tst, y_hat)^2,
    'mse' = mean((y_tst - y_hat)^2)
  )
  return(out)
}
phi_models_summary <- foreach(g = 1:ncol(phis), .combine = rbind) %do%
  phi_loop(g)


# Model would obviously do better if restricted to relevant Z's
# But is this consistent with the MLP architecture? And if so, how?







