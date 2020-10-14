# Set working directory
setwd('~/Documents/UCL')

# Load libraries, register cores
library(data.table)
library(ranger)
library(glmnet)
library(kernlab)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import, prep data
z <- readRDS('./grn/simulations/baseline.rds')
w <- readRDS('./grn/simulations/w.rds')
x <- readRDS('./grn/simulations/interventions.rds')
phis <- readRDS('./grn/simulations/phis.rds')
y <- readRDS('./grn/simulations/y.rds')

# Load lasso
trn_idx <- c(1:4000, 7001:1000)
tst_idx <- 4001:7000
phi_hat <- readRDS('./grn/phi_hat.rds')
lasso_f <- readRDS('./grn/lasso_f.rds')
y_hat <- predict(lasso_f, phi_hat[tst_idx, ], s = 'lambda.min')
lasso_df <- data.frame(
  expected = y[tst_idx], observed = as.numeric(y_hat)
)

################################################################################

### Benchmarks ###

# New W:
# For each of the following models:
# 1) lasso(Y|W,Z)
# 2) lasso(Y|Z,X)
# 3) SVM(Y|W,Z)
# 4) SVM(Y|Z,X)
# For each of {10, 20, ..., 100}% of the (training) data
# Test (20%) 3 models, report MSE.

n <- 3000
x1 <- cbind(w, z)
x2 <- cbind(z, x)
outer_loop <- function(b) {
  # Draw random test sample
  tst <- sample(tst_idx, 0.2 * n)
  # Permute training sample
  trn <- setdiff(tst_idx, tst)[sample.int(0.8 * n)]
  # Run inner loop
  inner_loop <- function(prop) {
    # Use first prop * length(trn) samples
    trn_p <- trn[seq_len(prop * length(trn))]
    # Fit E[Y|W,Z] models
    f1 <- cv.glmnet(x = x1[trn_p, ], y = y[trn_p])
    f2 <- ksvm(x = x1[trn_p, ], y = y[trn_p])
    # Fit E[Y|Z,X models]
    f3 <- cv.glmnet(x = x2[trn_p, ], y = y[trn_p])
    f4 <- ksvm(x = x2[trn_p, ], y = y[trn_p])
    # Test performance
    yhat_f1 <- as.numeric(predict(f1, x1[tst, ], s = 'lambda.min'))
    yhat_f2 <- as.numeric(predict(f2, x1[tst, ]))
    yhat_f3 <- as.numeric(predict(f3, x2[tst, ], s = 'lambda.min'))
    yhat_f4 <- as.numeric(predict(f4, x2[tst, ]))
    # Export
    out <- data.frame(
      'proportion' = prop, 'b' = b, 
      'model' = c('lasso_w,z', 'svm_w,z', 'lasso_z,x', 'svm_z,x', 'lasso_phi'),
      'mse' = c(
        mean((y[tst] - yhat_f1)^2), mean((y[tst] - yhat_f2)^2), 
        mean((y[tst] - yhat_f3)^2), mean((y[tst] - yhat_f4)^2), 
        with(lasso_df, mean((expected - observed)^2))
      )
    )
    return(out)
  }
  out <- foreach(p = seq(0.1, 1, 0.1), .combine = rbind) %dopar% 
    inner_loop(p)
  return(out)
}
# Execute in parallel
res <- foreach(i = 1:3, .combine = rbind) %do%
  outer_loop(i)

# Plot results
df <- as.data.table(res)
df[, err := mean(mse), by = .(proportion, model)]
df <- unique(select(df, -b, -mse))
colnames(df)[3] <- 'mse'
ggplot(df, aes(proportion, mse, color = model)) + 
  geom_line() + 
  scale_color_d3() +
  theme_bw()







# Which phis are TRUE
# Which are selected by lasso
# Which are selected by CI test




