#' Collapsed Gibbs sampler for infinite Bernoulli Mixture Model
#' 
#' Samples from the Dirichlet Process prior infinite BMM using
#' the Chinese Restaurant Process formulation
#' 
#' @inheritParams gibbs_full
#' @param a First Gamma parameter on prior on Alpha
#' @param b Second Gamma parameter on prior on Alpha
#' @param alpha Starting value of alpha
#' @param relabel Whether to apply Stephens 2000b relabelling method
#' @param burnrelabel If \code{relabel} is set then this parameter controls
#'   how many samples (from the end of the burnin period) are used for
#'   the initial batch run of Stephens' algorithm. Following this, the 
#'   online implementation is applied at each sample.
#' @param maxK Maximum number of clusters to assume. Reduces computational storage
#'   requirements.
#'   
#' @return A list with sampled components z, theta, alpha.
#'   If \code{relabel} is set then it also returns permutations,
#'   z_relabelled, and theta_relabelled
#'  
#' @export
gibbs_dp <- function(data, nsamples, alpha=NULL, a=1, b=1, beta=0.5, gamma=0.5, 
                     burnin=NULL, relabel=FALSE, burnrelabel=50, maxK=30, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    if (burnrelabel > burnin) burnrelabel <- round(0.1 * burnin)
    if (is.null(alpha)) alpha <- 0
    collapsed_gibbs_dp_cpp(data, nsamples, alpha, beta, gamma, a, b, 
                           burnin, relabel, burnrelabel, maxK, debug)
}

#' Collapsed Gibbs sampler for finite Bernoulli Mixture Model
#' 
#' @inheritParams gibbs_full
#' @return A list with sampled values z and theta.
#' @export
gibbs_collapsed <- function(data, nsamples, K, alpha=NULL, beta=0.5, gamma=0.5, 
                            a=1, b=1,
                            burnin=NULL, relabel=FALSE, burnrelabel=50, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    if (burnrelabel > burnin) burnrelabel <- round(0.1 * burnin)
    initial_K <- sample(1:K, nrow(data), replace=T)
    if (is.null(alpha)) alpha <- 0
    collapsed_gibbs_cpp(data, initial_K,
                        nsamples, K, alpha, beta, gamma, a, b, burnin, 
                        relabel, burnrelabel, debug)
}

#' Full Gibbs sampler for finite Bernoulli Mixture Model
#' 
#' @param data Data frame or matrix with observations in rows and binary
#'   variables in columns
#' @param nsamples Number of samples to take
#' @param K Number of mixtures
#' @param alpha Concentration parameter
#' @param beta First Beta parameter on prior used for all Bernoulli
#'   variables across all components
#' @param gamma Second Beta parameter on prior used for all Bernoulli
#'   variables across all components
#' @param burnin Number of samples to discard at start of chain
#' @param debug Whether to display debug messages.
#' @return A list with sampled values pi, z, and theta.
#' @export
gibbs_full <- function(data, nsamples, K, alpha=NULL, beta=0.5, gamma=0.5,
                       a=1, b=1,
                       burnin=NULL, relabel=FALSE, burnrelabel=50, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    initial_pi <- stats::runif(K)
    initial_pi <- exp(initial_pi)
    initial_pi <- initial_pi / sum(initial_pi)
    if (is.null(alpha)) alpha <- 0
    if (burnrelabel > burnin) burnrelabel <- round(0.1 * burnin)

    initial_theta <- matrix(stats::runif(K*ncol(data)), ncol=ncol(data), nrow=K)
    gibbs_cpp(data, initial_pi, initial_theta,
              nsamples, K, alpha, beta, gamma, a, b, burnin,
              relabel, burnrelabel, debug)
}

#' Stick breaking blocked Gibbs sampler for infinite Bernoulli Mixture Model
#' 
#' @param data Data frame or matrix with observations in rows and binary
#'   variables in columns
#' @param nsamples Number of samples to take
#' @param K Number of mixtures
#' @param alpha Concentration parameter
#' @param beta First Beta parameter on prior used for all Bernoulli
#'   variables across all components
#' @param gamma Second Beta parameter on prior used for all Bernoulli
#'   variables across all components
#' @param burnin Number of samples to discard at start of chain
#' @param debug Whether to display debug messages.
#' @return A list with sampled values pi, z, and theta.
#' @export
gibbs_stickbreaking <- function(data, nsamples, maxK, alpha=NULL, beta=0.5, gamma=0.5, a=1, b=1,
                       burnin=NULL, relabel=FALSE, burnrelabel=50, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    initial_pi <- stats::runif(maxK)
    initial_pi <- exp(initial_pi)
    initial_pi <- initial_pi / sum(initial_pi)
    if (is.null(alpha)) alpha <- 0

    initial_theta <- matrix(stats::runif(maxK*ncol(data)), ncol=ncol(data), nrow=maxK)
    gibbs_stickbreaking_cpp(data, initial_pi, initial_theta,
                            nsamples, maxK, alpha, beta, gamma, a, b, burnin, relabel, burnrelabel,
                            debug)
}

#' @export
alloc_test <- function() {
    library(BayesBinMix)
    set.seed(1)
    d <- 10 # number of columns
    n <- 50 # number of rows (sample size)
    K <- 2 	 # true number of clusters
    p.true <- myDirichlet(rep(10,K)) # true weight of each cluster
    z.true <- numeric(n) # true cluster membership
    z.true <- sample(K,n,replace=TRUE,prob = p.true)
    #true probability of positive responses per cluster:
    theta.true <- array(data = NA, dim = c(K,d)) 
    for(j in 1:d){
        theta.true[,j] <- rbeta(K, shape1 = 1, shape2 = 1)
    }
    x <- array(data=NA,dim = c(n,d)) # data: n X d array
    for(k in 1:K){
        myIndex <- which(z.true == k)
        for (j in 1:d){
            x[myIndex,j] <- rbinom(n = length(myIndex), 
                                   size = 1, prob = theta.true[k,j])   
        }
    }
    
    initial_K <- c(rep(1, nrow(x)/2),
                   rep(2, nrow(x)/2)) 
    foo <- allocation_cpp(x,  # df
                   initial_K,    # initial k
                   2,          # n samples
                   10,           # max K
                   1,            # alpha
                   0.5,          # beta
                   0.5,          # gamma
                   1,            # a
                   1,            # b
                   0,           # burnin
                   FALSE,        # relabel
                   1,           # burnrelabel
                   TRUE)         # debug
}

bbm_test <- function() {
    library(BayesBinMix)
    set.seed(1)
    d <- 10 # number of columns
    n <- 50 # number of rows (sample size)
    K <- 2 	 # true number of clusters
    p.true <- myDirichlet(rep(10,K)) # true weight of each cluster
    z.true <- numeric(n) # true cluster membership
    z.true <- sample(K,n,replace=TRUE,prob = p.true)
    #true probability of positive responses per cluster:
    theta.true <- array(data = NA, dim = c(K,d)) 
    for(j in 1:d){
        theta.true[,j] <- rbeta(K, shape1 = 1, shape2 = 1)
    }
    x <- array(data=NA,dim = c(n,d)) # data: n X d array
    for(k in 1:K){
        myIndex <- which(z.true == k)
        for (j in 1:d){
            x[myIndex,j] <- rbinom(n = length(myIndex), 
                                   size = 1, prob = theta.true[k,j])   
        }
    }
    foo <- allocationSamplerBinMix2(Kmax=10, m=1000, burn=10, data=x,
                                    outputDir = "foo", ClusterPrior = "poisson",
                                    zStart = c(rep(1, nrow(x)/2), rep(2, nrow(x)/2)))
}


#' Plots samples from Bernoulli mixture models
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Rcpp evalCpp
#' @export
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
            ggplot2::ggplot(ggplot2::aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                ggplot2::geom_line() +
                ggplot2::theme_bw() +
                ggplot2::labs(x="Sample", y="Pi") +
                ggplot2::scale_colour_discrete("Cluster")
        plts[[length(plts) + 1]] <- plt_pi
    }

    if (z) {
        if (is.null(cluster_labels)) {
            cluster_labels <- 1:K
        }
        dimnames(z_raw) <- list('sample'=1:S, 'observation'=1:N)
        z_long <- as.data.frame.table(z_raw, responseName="cluster") %>%
                    dplyr::mutate(cluster = factor(cluster, levels=1:K, labels=cluster_labels))

        z_props <- z_long %>%
            dplyr::filter(sample != 1) %>%
            dplyr::group_by(sample, cluster) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::mutate(prop = n / sum(n))

        cluster_to_plot <- z_props %>%
            dplyr::filter(prop > cluster_threshold) %>%
            dplyr::distinct(sample, cluster)
        unique_clusters <- unique(cluster_to_plot$cluster)

        plt_z <- z_props %>%
            dplyr::mutate(cluster = factor(cluster, levels=unique_clusters)) %>%
            ggplot2::ggplot(ggplot2::aes(x=as.integer(sample), y=prop, colour=cluster)) +
                ggplot2::geom_line() +
                ggplot2::theme_bw() +
                ggplot2::ylim(0, 1) +
                ggplot2::labs(x="Sample", y="Proportion in cluster") +
                ggplot2::scale_colour_discrete("Cluster", guide=F, drop=F)
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
                        dplyr::mutate(cluster = factor(cluster, levels=cluster_labels),
                               theta_var = factor(theta_var, levels=theta_labels))

        if (!is.null(theta_to_display)) {
            theta_long <- theta_long %>%
                        dplyr::filter(theta_var %in% theta_to_display) %>%
                        dplyr::mutate(theta_var = factor(theta_var, levels=theta_to_display))
        }

        foo <- cluster_to_plot %>%
            dplyr::left_join(theta_long, by=c('cluster'='cluster', 'sample'='sample'))
        plt_theta <- foo %>%
            dplyr::filter(sample != 1) %>%
            dplyr::mutate(cluster = factor(cluster, levels=unique_clusters)) %>%
            ggplot2::ggplot(ggplot2::aes(x=as.integer(sample), y=value, colour=as.factor(cluster))) +
                ggplot2::geom_line() +
                ggplot2::facet_wrap(~theta_var) +
                ggplot2::theme_bw() +
                ggplot2::ylim(0, 1) +
                ggplot2::labs(x="Sample", y="Theta") +
                ggplot2::scale_colour_discrete("Cluster", guide=F, drop=F)
        plts[[length(plts) + 1]] <- plt_theta
    }
    gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs=plts, ncol=1, heights=heights))
}

plot_alpha <- function(obj) {
    ggplot2::ggplot(data.frame(foo=obj$alpha), ggplot2::aes(foo)) +
        ggplot2::geom_histogram(binwidth = 0.1, colour="black", fill="white") +
        ggplot2::xlim(0, 5)
}
