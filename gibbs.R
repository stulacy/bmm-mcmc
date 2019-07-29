# Runs Gibbs sampling
# Starting with fixed K
# Then will buildup to random K using Dirichlet prior
# Then will finalise with infinite mixture using Dirichlet Process prior
library(tidyverse)
library(gridExtra)
source("utils.R")

# Ok this has seemed to work on an easy dataset with 100 observations
# and 2 well separated classes
df <- readRDS("data/K2_N100_P5_clean.rds")

# Form dataset that know what the first values should be
# N = 4, P = 3, K =2
test_df <- matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0), nrow=4, ncol=3, byrow=T)
test_df

# debug
set.seed(12)
samples_cpp <- gibbs_collapsed(df[1:7, ], 2, K=2, debug=TRUE)

set.seed(12)
foo <- gibbs_dp(df[1:7, ], 5, debug=TRUE)

# What about on the same dataset with 1 thousand observations?
df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
samples <- gibbs_collapsed(df_2, 1000, K=2)
plot_gibbs(samples)

samples_dp <- gibbs_dp(df_2, 10000, debug=FALSE, alpha = 0.1)
plot_gibbs(samples_dp, pi=F, cluster_threshold = 0.1)

# Ok so seems to be fine with the number of observations, indeed it found
# N=1000 much easier than N=100

# So is it the number of clusters that's the problem?
# Let's try using K=3
# Oh it does seem to have worked now have separated clusters more
df_3 <- readRDS("data/K3_N1000_P5_clean.rds")
samples <- gibbs_collapsed(df_3, 1000, K=3)
plot_gibbs(samples)

samples_dp3 <- gibbs_dp(df_3, 20000)
plot_gibbs(samples_dp3, cluster_threshold = 0.10)

# Testing full Gibbs sampling and can see that like with the Collapsed Gibbs,
# it works fine in the situation with K=2, N=1000.
# And furthermore can easily obtain thetas, which must be obtainable from
# collapsed gibbs sampler but I just don't know how.
foo <- gibbs_full(df_2, 1000, 2, debug=FALSE)
plot_gibbs(foo)

# Can it handle K=3 however?
# Yes it can rather easily
foo <- gibbs_full(df_3, 1000, 3)
plot_gibbs(foo)
