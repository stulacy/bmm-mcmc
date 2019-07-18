# Runs Gibbs sampling
# Starting with fixed K
# Then will buildup to random K using Dirichlet prior
# Then will finalise with infinite mixture using Dirichlet Process prior
library(tidyverse)
library(gridExtra)

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
    list(z=allocations)
}

gibbs_collapsed_cpp_wrapper <- function(df, nsamples, K, alpha=1, beta=0.5, gamma=0.5, debug=FALSE) {
    initial_K <- t(rmultinom(nrow(df), 1, rep(1/K, K)))
    collapsed_gibbs_cpp(df, initial_K,
                        nsamples, K, alpha, beta, gamma, debug)
}

gibbs_full_cpp_wrapper <- function(data, nsamples, K, alpha=1, beta=0.5, gamma=0.5,
                                   debug=FALSE) {
    initial_pi <- runif(K)
    initial_pi <- exp(initial_pi)
    initial_pi <- initial_pi / sum(initial_pi)

    initial_theta <- matrix(runif(K*ncol(data)), ncol=ncol(data), nrow=K)
    gibbs_cpp(data, initial_pi, initial_theta,
              nsamples, K, alpha, beta, gamma, debug)
}

plot_gibbs_collapsed <- function(obj) {
    z <- obj$z
    N <- dim(z)[1]
    K <- dim(z)[2]
    S <- dim(z)[3]
    dimnames(z) <- list('observation'=1:N, 'cluster'=1:K, 'sample'=1:S)
    z_long <- as.data.frame.table(z, responseName="value")
    
    z_long %>% 
        filter(sample != 1) %>%
        group_by(sample, cluster) %>%
        summarise(n = sum(value)) %>%
        mutate(prop = n / sum(n)) %>%
        ggplot(aes(x=as.integer(sample), y=prop, colour=as.factor(cluster))) +
            geom_line(alpha=0.5) +
            theme_bw() +
            labs(x="Sample", y="Proportion in cluster") +
            scale_colour_discrete("Cluster")
}

plot_gibbs_complete <- function(obj, theta=TRUE, z=TRUE, pi=TRUE, heights=NULL) {

    plts <- list()

    if (pi) {
        pi <- obj$pi
        S <- dim(pi)[1]
        K <- dim(pi)[2]
        dimnames(pi) <- list('sample'=1:S, 'cluster'=1:K)
        pi_long <- as.data.frame.table(pi, responseName="value")

        plt_pi <- pi_long %>%
            ggplot(aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                geom_line(alpha=0.5) +
                theme_bw() +
                labs(x="Sample", y="Pi") +
                scale_colour_discrete("Cluster")
        plts[[length(plts) + 1]] <- plt_pi
    }

    if (z) {
        z <- obj$z
        N <- dim(z)[1]
        K <- dim(z)[2]
        S <- dim(z)[3]
        dimnames(z) <- list('observation'=1:N, 'cluster'=1:K, 'sample'=1:S)
        z_long <- as.data.frame.table(z, responseName="value")

        plt_z <- z_long %>%
            filter(sample != 1) %>%
            group_by(sample, cluster) %>%
            summarise(n = sum(value)) %>%
            mutate(prop = n / sum(n)) %>%
            ggplot(aes(x=as.integer(sample), y=prop, colour=as.factor(cluster))) +
                geom_line(alpha=0.5) +
                theme_bw() +
                labs(x="Sample", y="Proportion in cluster") +
                scale_colour_discrete("Cluster")
        plts[[length(plts) + 1]] <- plt_z
    }

    if (theta) {
        theta <- obj$theta
        K <- dim(theta)[1]
        P <- dim(theta)[2]
        S <- dim(theta)[3]
        dimnames(theta) <- list('cluster'=1:K, 'variable'=1:P, 'sample'=1:S)
        theta_long <- as.data.frame.table(theta, responseName = "value")

        plt_theta <- theta_long %>%
            ggplot(aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                geom_line(alpha=0.5) +
                facet_wrap(~variable) +
                theme_bw() +
                labs(x="Sample", y="Theta") +
                scale_colour_discrete("Cluster")
        plts[[length(plts) + 1]] <- plt_theta
    }
    grid.arrange(arrangeGrob(grobs=plts, ncol=1, heights=heights))
}

# Ok this has seemed to work on an easy dataset with 100 observations
# and 2 well separated classes
df <- readRDS("data/K2_N100_P5_clean.rds")
samples_R <- gibbs_collapsed(df, 100, K=2)
plot_gibbs_collapsed(samples_R)

set.seed(12)
samples_cpp <- gibbs_collapsed_cpp_wrapper(df[1:5, ], 3, K=2, debug = TRUE)
plot_gibbs_collapsed(samples_cpp)
samples_cpp$theta

# What about on the same dataset with 1 thousand observations?
df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_2, 1000, K=2)
plot_gibbs_collapsed(samples)

# Ok so seems to be fine with the number of observations, indeed it found
# N=1000 much easier than N=100

# So is it the number of clusters that's the problem?
# Let's try using K=3
# Oh it does seem to have worked now have separated clusters more
df_3 <- readRDS("data/K3_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_3, 10000, K=3)
plot_gibbs_collapsed(samples)

# Testing full Gibbs sampling and can see that like with the Collapsed Gibbs,
# it works fine in the situation with K=2, N=1000.
# And furthermore can easily obtain thetas, which must be obtainable from
# collapsed gibbs sampler but I just don't know how.
foo <- gibbs_full_cpp_wrapper(df_2, 1000, 2, debug=FALSE)
plot_gibbs_complete(foo)

# Can it handle K=3 however?
# Yes it can rather easily
foo <- gibbs_full_cpp_wrapper(df_3, 10000, 3)
plot_gibbs_complete(foo)
