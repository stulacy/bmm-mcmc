## Runs Gibbs sampling
## Starting with fixed K
## Then will buildup to random K using Dirichlet prior
## Then will finalise with infinite mixture using Dirichlet Process prior
#library(tidyverse)
#library(gridExtra)
#source("R/utils.R")
#
#df <- readRDS("data/K2_N100_P5_clean.rds")
#df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
#df_3 <- readRDS("data/K3_N1000_P5_clean.rds")
#
## Form dataset that know what the first values should be
## N = 4, P = 3, K =2
#test_df <- matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0), nrow=4, ncol=3, byrow=T)
#test_df
#
## debug
#set.seed(12)
#samples_cpp <- gibbs_collapsed(df[1:7, ], 2, K=2, debug=TRUE)
#
#set.seed(12)
#foo <- gibbs_dp(df[1:7, ], 5, debug=TRUE)
#
## What about on the same dataset with 1 thousand observations?
#samples <- gibbs_collapsed(df_2, 1000, K=2)
#plot_gibbs(samples)
#
#samples_dp <- gibbs_dp(df_2, 10000, debug=FALSE, alpha = 0.1)
#plot_gibbs(samples_dp, cluster_threshold = 0.1)
#
## Ok so seems to be fine with the number of observations, indeed it found
## N=1000 much easier than N=100
#
## So is it the number of clusters that's the problem?
## Let's try using K=3
## Oh it does seem to have worked now have separated clusters more
#samples <- gibbs_collapsed(df_3, 1000, K=3)
#plot_gibbs(samples)
#
#samples_dp3 <- gibbs_dp(df_3, 20000)
#plot_gibbs(samples_dp3, cluster_threshold = 0.10)
#
## Testing full Gibbs sampling and can see that like with the Collapsed Gibbs,
## it works fine in the situation with K=2, N=1000.
## And furthermore can easily obtain thetas, which must be obtainable from
## collapsed gibbs sampler but I just don't know how.
#foo <- gibbs_full(df_2, 1000, 2, debug=FALSE)
#plot_gibbs(foo)
#
## Can it handle K=3 however?
## Yes it can rather easily
#foo <- gibbs_full(df_3, 1000, 3)
#plot_gibbs(foo)
#
#
######### Label switching
#samples_label <- gibbs_collapsed(df_3, K=20, 1000, burnin=500)
## Annoyingly I can't force data to have label switching problem
## But will go through the code anyway, then can test on DP model which definitely
## exhibits label switching
#plot_gibbs(samples_label)
#
#library(label.switching)
#
## Testing all 3 implementations. One is Stephens 2000
#labs_steph <- stephens(samples_label$probabilities)
## Other 2 are ECR from Papastamoulis 2010
## First just requires labels
#labs_ecr1 <- ecr.iterative.1(samples_label$z, dim(samples_label$theta)[1], opt_init = 1:20)
## Second requires probabilities
#labs_ecr2 <- ecr.iterative.2(samples_label$z, dim(samples_label$probabilities)[3], samples_label$probabilities)
#
## They don't fully agree, although that could be because I'm using a silly example here forcing k=10
#all(labs_steph$permutations == labs_ecr1$permutations)
#all(labs_steph$permutations == labs_ecr2$permutations)
#all(labs_ecr1$permutations == labs_ecr2$permutations)
#
## How to use these
#reorder_clusters <- function(obj, perms) {
#    nsamps <- dim(obj$z)[1]
#    D <- dim(obj$theta)[2]
#    for (j in 1:nsamps) {
#        for (d in 1:D) {
#            obj$theta[, d, j] <- obj$theta[perms[j, ], d, j]
#        }
#    }
#    if (!is.null(obj$probabilities)) {
#        N <- dim(obj$probabilities)[2]
#        for (j in 1:nsamps) {
#            for (i in 1:N) {
#                obj$probabilities[j, i, ] <- obj$probabilities[j, i, perms[j, ]]
#            }
#        }
#    }
#    obj
#}
#
#plot_gibbs(samples_label)
#plot_gibbs(reorder_clusters(samples_label, labs_steph$permutations))
## Looks like ECR1 has failed,i.e. need probabilities
#plot_gibbs(reorder_clusters(samples_label, labs_ecr1$permutations))
#plot_gibbs(reorder_clusters(samples_label, labs_ecr2$permutations))
#
## However the iterative ECR method is far quicker so let's stick with that
#microbenchmark::microbenchmark(stephens(samples_label$probabilities),
#                               ecr.iterative.1(samples_label$z, dim(samples_label$theta)[1]),
#                               ecr.iterative.2(samples_label$z, 
#                                               dim(samples_label$probabilities)[3],
#                                               samples_label$probabilities),
#                               times=3)
#
##### Now since I'm not seeing the relabelling issue on finite K
## let's see if we can fix the Dirichlet Process Model where it does occur
## NB: Since I can't currently save the cluster probabilities for this model
## I'm going to use the ECR1 method instead, which just uses hard labels
#samples_dp <- gibbs_dp(df_3, 1000, burnin=500)
## Not a huge amount of label switching but still definitely present
#plot_gibbs(samples_dp)
#
#labs_dp <- ecr.iterative.1(samples_dp$z, dim(samples_dp$theta)[1])
#
## They don't fully agree, although that could be because I'm using a silly example here forcing k=10
#all(labs_steph$permutations == labs_ecr1$permutations)
#all(labs_steph$permutations == labs_ecr2$permutations)
#all(labs_ecr1$permutations == labs_ecr2$permutations)
#
## How to use these
#reorder_clusters <- function(obj, perms) {
#    nsamps <- dim(obj$z)[1]
#    D <- dim(obj$theta)[2]
#    for (j in 1:nsamps) {
#        for (d in 1:D) {
#            obj$theta[, d, j] <- obj$theta[perms[j, ], d, j]
#        }
#    }
#    if (!is.null(obj$probabilities)) {
#        N <- dim(obj$probabilities)[2]
#        for (j in 1:nsamps) {
#            for (i in 1:N) {
#                obj$probabilities[j, i, ] <- obj$probabilities[j, i, perms[j, ]]
#            }
#        }
#    }
#    obj
#}
#
#plot_gibbs(samples_label)
#plot_gibbs(reorder_clusters(samples_label, labs_steph$permutations))
## Looks like ECR1 has failed,i.e. need probabilities
#plot_gibbs(reorder_clusters(samples_label, labs_ecr1$permutations))
#plot_gibbs(reorder_clusters(samples_label, labs_ecr2$permutations))
#