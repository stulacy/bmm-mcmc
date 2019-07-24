# Runs Gibbs sampling
# Starting with fixed K
# Then will buildup to random K using Dirichlet prior
# Then will finalise with infinite mixture using Dirichlet Process prior
library(tidyverse)
library(gridExtra)

gibbs_dp_cpp_wrapper <- function(df, nsamples, alpha=1, beta=0.5, gamma=0.5, debug=FALSE) {
    collapsed_gibbs_dp_cpp(df, nsamples, alpha, beta, gamma, debug)
}

gibbs_collapsed_cpp_wrapper <- function(df, nsamples, K, alpha=1, beta=0.5, gamma=0.5, debug=FALSE) {
    initial_K <- sample(1:K, nrow(df), replace=T)-1
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

plot_gibbs <- function(obj, theta=TRUE, z=TRUE, pi=FALSE, heights=NULL, cluster_threshold=0.1,
                       cluster_labels=NULL, theta_labels=NULL, theta_to_display=NULL) {

    plts <- list()

    theta_raw <- obj$theta
    z_raw <- obj$z

    K <- dim(theta_raw)[1]
    P <- dim(theta_raw)[2]
    S <- dim(theta_raw)[3]
    N <- ncol(z_raw)

    if (pi) {
        pi <- obj$pi
        S <- dim(pi)[1]
        K <- dim(pi)[2]
        if (is.null(cluster_labels)) {
            cluster_labels <- 1:K
        }
        dimnames(pi) <- list('sample'=1:S, 'cluster'=cluster_labels)
        pi_long <- as.data.frame.table(pi, responseName="value")

        plt_pi <- pi_long %>%
            ggplot(aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                geom_line() +
                theme_bw() +
                labs(x="Sample", y="Pi") +
                scale_colour_discrete("Cluster")
        plts[[length(plts) + 1]] <- plt_pi
    }

    if (z) {
        if (is.null(cluster_labels)) {
            cluster_labels <- 1:K
        }
        dimnames(z_raw) <- list('sample'=1:S, 'observation'=1:N)
        z_long <- as.data.frame.table(z_raw, responseName="cluster") %>%
                    mutate(cluster = factor(cluster, levels=1:K, labels=cluster_labels))

        z_props <- z_long %>%
            filter(sample != 1) %>%
            group_by(sample, cluster) %>%
            summarise(n = n()) %>%
            mutate(prop = n / sum(n))

        cluster_to_plot <- z_props %>%
            filter(prop > cluster_threshold) %>%
            distinct(sample, cluster)

        plt_z <- z_props %>%
            ggplot(aes(x=as.integer(sample), y=prop, colour=cluster)) +
                geom_line() +
                theme_bw() +
                ylim(0, 1) +
                labs(x="Sample", y="Proportion in cluster") +
                scale_colour_discrete("Cluster", guide=F, drop=F)
        plts[[length(plts) + 1]] <- plt_z
    }

    if (theta) {
        if (is.null(cluster_labels)) {
            cluster_labels <- 1:K
        }
        if (is.null(theta_labels)) {
            theta_labels <- 1:P
        }
        dimnames(theta_raw) <- list('cluster'=cluster_labels, 'theta_var'=theta_labels, 'sample'=1:S)
        theta_long <- as.data.frame.table(theta_raw, responseName = "value") %>%
                        mutate(cluster = factor(cluster, levels=cluster_labels),
                               theta_var = factor(theta_var, levels=theta_labels))

        if (!is.null(theta_to_display)) {
            theta_long <- theta_long %>%
                            filter(theta_var %in% theta_to_display)
        }

        foo <- cluster_to_plot %>%
            left_join(theta_long, by=c('cluster'='cluster', 'sample'='sample'))
        plt_theta <- foo %>%
            filter(sample != 1) %>%
            ggplot(aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                geom_line() +
                facet_wrap(~theta_var) +
                theme_bw() +
                ylim(0, 1) +
                labs(x="Sample", y="Theta") +
                scale_colour_discrete("Cluster", guide=F, drop=F)
        plts[[length(plts) + 1]] <- plt_theta
    }
    grid.arrange(arrangeGrob(grobs=plts, ncol=1, heights=heights))
}

# Ok this has seemed to work on an easy dataset with 100 observations
# and 2 well separated classes
df <- readRDS("data/K2_N100_P5_clean.rds")
#samples_R <- gibbs_collapsed(df, 100, K=2)
#plot_gibbs(samples_R)

# Form dataset that know what the first values should be
# N = 4, P = 3, K =2
test_df <- matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0), nrow=4, ncol=3, byrow=T)
test_df
initial_K <- matrix(c(1, 0, 0, 1, 1, 0, 0, 1), nrow=4, ncol=2, byrow=T)
initial_K

collapsed_gibbs_cpp(test_df, initial_K, 2, 2, 1, 0.5, 0.5, TRUE)

# debug
set.seed(12)
gibbs_dp_cpp_wrapper(df[1:7, ], 2, debug=TRUE)

set.seed(12)
samples_cpp <- gibbs_collapsed_cpp_wrapper(df[1:7, ], 2, K=2, debug=TRUE)

plot_gibbs(samples_cpp, pi=F)

# What about on the same dataset with 1 thousand observations?
df_2 <- readRDS("data/K2_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_2, 1000, K=2)
plot_gibbs(samples, pi=F))

samples_dp <- gibbs_dp_cpp_wrapper(df_2, 10000, debug=FALSE)
plot_gibbs(samples_dp, pi=F, cluster_threshold = 0.1)

# Ok so seems to be fine with the number of observations, indeed it found
# N=1000 much easier than N=100

# So is it the number of clusters that's the problem?
# Let's try using K=3
# Oh it does seem to have worked now have separated clusters more
df_3 <- readRDS("data/K3_N1000_P5_clean.rds")
samples <- gibbs_collapsed_cpp_wrapper(df_3, 1000, K=3)
plot_gibbs(samples, pi=F)

samples_dp3 <- gibbs_dp_cpp_wrapper(df_3, 1000)
plot_gibbs_dp(samples_dp3, cluster_threshold = 0.15)

# Testing full Gibbs sampling and can see that like with the Collapsed Gibbs,
# it works fine in the situation with K=2, N=1000.
# And furthermore can easily obtain thetas, which must be obtainable from
# collapsed gibbs sampler but I just don't know how.
foo <- gibbs_full_cpp_wrapper(df_2, 1000, 2, debug=FALSE)
plot_gibbs(foo)

# Can it handle K=3 however?
# Yes it can rather easily
foo <- gibbs_full_cpp_wrapper(df_3, 1000, 3)
plot_gibbs(foo)
