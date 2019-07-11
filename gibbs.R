# Runs Gibbs sampling
# Starting with fixed K
# Then will buildup to random K using Dirichlet prior
# Then will finalise with infinite mixture using Dirichlet Process prior
library(tidyverse)

# Collapsed Gibbs sampler
# For each i in 1:N
#  Sample z_i given (theta, z_(-i), x)
# For k in 1:K
#  Sample theta_k given (theta_(-k), z, x)

# Tu in Section 3.23 says needs to sample p(z | ...), p(pi | ...), and p(theta | ...)
# But if use collapsed Gibbs then can sample p(z) directly
# https://people.eecs.berkeley.edu/~stephentu/writeups/mixturemodels.pdf

# Eq 21 from Van Maarten also suggests can just use collapsed Gibbs sampler over z
# this equation looks to be the same as Tu's version.
# https://pdfs.semanticscholar.org/525c/ff658d34ae4d47e84b8ec4ede3ce6c561afc.pdf

# And this looks to be the same as 3.5 from Neal
# http://www.stat.columbia.edu/npbayes/papers/neal_sampling.pdf
gibbs_collapsed <- function(df, nsamples, K, alpha=1, beta=0.5, gamma=0.5, verbose=FALSE) {
    
    N <- nrow(df)
    P <- ncol(df)

    # Random sample for first row
    allocations <- matrix(0, nrow=nsamples, ncol=N)
    allocations[1, ] <- sample(1:K, N, replace=T)
    for (j in 2:nsamples) {
        cat(sprintf("\nSample: %d\n ", j))
        # Set params
        for (i in 1:N) {
            probs <- rep(0, K)
            for (k in 1:K) {
                # All points in current cluster except current data point. HOW DOES THIS DIFFER FROM N_nk?
                ck <- allocations[j-1, -i] == k
                N_nk <- sum(ck)
    
                frac <- log(N_nk + alpha/K) - log(N - 1 + alpha)
                if (verbose) cat(sprintf("Frac: %f\n", frac))
    
                loglh <- 0
                for (d in 1:P) {
                    sum_x <- sum(df[-i, ][ck, d])
    
                    num_1 <- df[i,d] * log(beta + sum_x)
                    num_2 <- (1-df[i,d]) * log(gamma + N_nk - sum_x)
                    denom <- log(beta + gamma + N_nk)
                    contrib <- num_1 + num_2 - denom
                    loglh <- loglh + contrib
                }
                if (verbose) cat(sprintf("loglh: %f\n", loglh))
                probs[k] <- frac + loglh
            }
            if (verbose) cat(sprintf("raw probs: %f\n", probs))
            probs <- exp(probs)
            sum_ <- sum(probs)
            probs <- probs / sum_
            if (verbose) cat(sprintf("normalised probs: %f\n", probs))
            #cat(paste0("Probs: ", probs, "\n"))
            allocations[j, i] <- sample(1:k, 1, prob=probs)
        }
    }
    allocations
}

gibbs_collapsed_cpp_wrapper <- function(df, nsamples, K, alpha=1, beta=0.5, gamma=0.5, verbose=FALSE) {
    initial_K <- sample(1:K, nrow(df), replace=T)
    collapsed_gibbs_cpp(df, initial_K,
                        nsamples, K, alpha, beta, gamma, verbose)
}
    

plot_gibbs <- function(samples) {
    K <- length(unique(samples[1, ]))
    nsamples <- nrow(samples)
    nobs <- ncol(samples)
    foo <- data.frame(t(apply(samples, 1, function(row) {
                              sapply(1:K, function(x) sum(row == x))       
    })))
    colnames(foo) <- paste0("Cluster", seq(K))
    foo$step <- 1:nsamples
    
    foo %>%
        gather(cluster, num, -step) %>%
        mutate(prop = num/nobs) %>%
        ggplot(aes(x=step, y=prop, colour=cluster)) +
            geom_line(alpha=0.5) +
            scale_colour_discrete("") +
            labs(x="Step number", y="Proportion") +
            theme_bw()
}

# Ok this has seemed to work on an easy dataset with 100 observations 
# and 2 well separated classes
df <- readRDS("data/K2_N100_P5_clean.rds")
samples <- gibbs_collapsed(df, 1000, K=2)
plot_gibbs(samples)
samples_cpp <- gibbs_collapsed_cpp_wrapper(df, 1000, K=2, verbose = FALSE)
plot_gibbs(samples_cpp)

# What about on the same dataset with 1 thousand observations?
df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_2, 1000, K=2)
plot_gibbs(samples)

# Ok so seems to be fine with the number of observations, indeed it found
# N=1000 much easier than N=100

# So is it the number of clusters that's the problem?
# Let's try using K=3
# Yep seems to be affected by the label switching problem
df_3 <- readRDS("data/K3_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_3, 10000, K=3)
plot_gibbs(samples)

##### TODO:
#  - Get working with multiple clusters
#  - Extract Bernoulli parameters
#  - Try getting it working with redundant features
#  - Dirichlet process