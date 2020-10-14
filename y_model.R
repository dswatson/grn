# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(ranger)
library(glmnet)
library(kernlab)
library(broom)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import, prep data
z <- readRDS('./grn/simulations/baseline.rds')
w <- readRDS('./grn/simulations/w.rds')
x <- readRDS('./grn/simulations/interventions.rds')
phis <- readRDS('./grn/simulations/phis.rds')
k <- ncol(phis)

# Most/least variable in W?
rho_fn <- function(phi_idx) {
  cor.test(phis[, phi_idx], w, method = 'spearman') %>%
    tidy(.) %>%
    mutate(idx = phi_idx) %>%
    select(idx, statistic, p.value)
}
df <- foreach(j = seq_len(k), .combine = rbind) %dopar%
  rho_fn(j) %>%
  mutate(q.value = p.adjust(p.value, method = 'fdr')) %>%
  arrange(p.value)

# Desiderata for Y: 
# (1) linear combination of phis such that
# (2) some phis receive zero weight,
# (3) some phis with nonzero weight are invariant w/r/t W, and
# (4) some phis with nonzero weight vary w/r/t W -- i.e., are causal mediators.

# Define beta vector
most_variable <- head(df$idx, 25)
least_variable <- tail(df$idx, 25)
amplitude <- 4
nonzeros <- rnorm(50, mean = amplitude, sd = 1) * 
  sample(c(1, -1), size = 50, replace = TRUE)
beta <- double(k)
beta[c(most_variable, least_variable)] <- nonzeros

# Define response
n <- nrow(phis)
y <- phis %*% beta + rnorm(n, sd = 1)
saveRDS(y, './grn/simulations/y.rds')

################################################################################

### LEARN LASSO WEIGHTS ###

# Import adjacency matrix
imp <- fread('./grn/adj_mat.csv')
adj_mat <- ifelse(imp >= 10, 1, 0)
outdegree <- colSums(adj_mat)
keep <- which(outdegree >= 100)

# Propagate phis
phi_predictor <- function(phi_idx) {
  # Prep data
  tf <- keep[phi_idx]
  x <- cbind(z[, c(tf, 334 + which(adj_mat[, tf] == 1))], w)
  # Import model, export phi
  f <- readRDS(paste0('./grn/phi_cont/phi', phi_idx, '.rds'))
  phi_j_hat <- predict(f, x)$predictions
  return(phi_j_hat)
}
phi_hat <- foreach(j = seq_len(k), .combine = cbind) %dopar% 
  phi_predictor(j)
colnames(phi_hat) <- paste0('phi', seq_len(k))
saveRDS(phi_hat, './grn/phi_hat.rds')

# Learn lasso weights
trn_idx <- c(1:4000, 7001:1000)
tst_idx <- 4001:7000
lasso_f <- cv.glmnet(x = phi_hat[trn_idx, ], y = y[trn_idx], parallel = TRUE)

# Coefficients look ok?
beta_df <- data.frame(
  expected = beta, observed = as.numeric(coef(lasso_f, s = 'lambda.min'))[-1]
)
lo <- min(beta_df$expected, beta_df$observed)
hi <- max(beta_df$expected, beta_df$observed)
ggplot(beta_df, aes(expected, observed)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  xlim(lo, hi) + ylim(lo, hi) + 
  theme_bw()
with(beta_df, cor(expected, observed)^2)

# Model fit?
y_hat <- predict(lasso_f, phi_hat[tst_idx, ], s = 'lambda.min')
lasso_df <- data.frame(
  expected = y[tst_idx], observed = as.numeric(y_hat)
)
lo <- min(lasso_df$expected, lasso_df$observed)
hi <- max(lasso_df$expected, lasso_df$observed)
ggplot(lasso_df, aes(expected, observed)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  xlim(lo, hi) + ylim(lo, hi) + 
  theme_bw()
with(lasso_df, cor(expected, observed)^2)
saveRDS(lasso_f, './grn/lasso_f.rds')







