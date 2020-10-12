# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(ranger)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(42, kind = "L'Ecuyer-CMRG")

# Import, prep data
x_dat <- readRDS('./grn/simulations/baseline.rds') %>%
  as.data.table(.) %>%
  mutate(W = paste0('W', rep(1:10, each = 2000)))
phis <- readRDS('./grn/simulations/phis.rds')

# Random forest loop
phi_rf_loop <- function(tf) {
  idx <- sample.int(20000, 2000)
  f <- ranger(x = x_dat[idx, ], y = phis[idx, tf], mtry = 1000,
              num.trees = 500, num.threads = 8)
  saveRDS(f, paste0('./grn/phi_models/phi', tf, '.rds'))
  out <- data.table(
    'phi' = paste0('phi', tf),
     'r2' = f$r.squared,
    'mse' = f$prediction.error
  )
  return(out)
}
phi_models_summary <- foreach(g = 1:ncol(phis), .combine = rbind) %do%
  phi_rf_loop(g)


# Model would obviously do better if restricted to relevant Z's
# But is this consistent with the MLP architecture? And if so, how?







