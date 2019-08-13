################################# Compile C++ source and load utility functions
library(tidyverse)
library(gridExtra)
source("utils.R")

# Load simulated dataset
df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
# This has produced 2 clusters with ratio 0.7 : 0.3 with 5 variables with following thetas
# [0.7, 0.8, 0.2, 0.1, 0.1,   # cluster 1 theta
#  0.2, 0.2, 0.9, 0.8, 0.6]   # cluster 2 theta

n_samples <- 10000

# Full Gibbs sampler with Finite K
samples_full <- gibbs_full(df_2, n_samples, K=10)
plot_gibbs(samples_full)

# Collapsed Gibbs sampler with Finite K
samples_collapsed <- gibbs_collapsed(df_2, n_samples, K=2)
plot_gibbs(samples_collapsed)

# Collapsed Gibbs Dirichlet Process with infinite K
samples_dp <- gibbs_dp(df_2, n_samples, debug=FALSE)
plot_gibbs(samples_dp, cluster_threshold = 0.1)

# Collapsed Gibbs sampler with Finite K=10
samples_collapsed10 <- gibbs_collapsed(df_2, n_samples, K=10)
plot_gibbs(samples_collapsed10)
