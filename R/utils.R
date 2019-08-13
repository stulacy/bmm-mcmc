plot_alpha <- function(obj) {
    ggplot2::ggplot(data.frame(foo=obj$alpha), ggplot2::aes(foo)) +
        ggplot2::geom_histogram(binwidth = 0.1, colour="black", fill="white") +
        ggplot2::xlim(0, 5)
}

gibbs_dp <- function(df, nsamples, a=1, b=1, alpha=1, beta=0.5, gamma=0.5, 
                     burnin=NULL, relabel=TRUE, burnrelabel=50, maxK=30, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    if (burnrelabel > burnin) burnrelabel <- round(0.1 * burnin)
    collapsed_gibbs_dp_cpp(df, nsamples, alpha, beta, gamma, a, b, 
                           burnin, relabel, burnrelabel, maxK, debug)
}

gibbs_collapsed <- function(df, nsamples, K, alpha=1, beta=0.5, gamma=0.5, burnin=NULL, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    initial_K <- sample(1:K, nrow(df), replace=T)
    collapsed_gibbs_cpp(df, initial_K,
                        nsamples, K, alpha, beta, gamma, burnin, debug)
}

gibbs_full <- function(data, nsamples, K, alpha=1, beta=0.5, gamma=0.5,
                       burnin=NULL, debug=FALSE) {
    if (is.null(burnin)) burnin <- round(0.1 * nsamples)
    initial_pi <- runif(K)
    initial_pi <- exp(initial_pi)
    initial_pi <- initial_pi / sum(initial_pi)

    initial_theta <- matrix(runif(K*ncol(data)), ncol=ncol(data), nrow=K)
    gibbs_cpp(data, initial_pi, initial_theta,
              nsamples, K, alpha, beta, gamma, burnin, debug)
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
            gglot2::ggplot(gglot2::aes(x=as.integer(sample), y=prop, colour=cluster)) +
                gglot2::geom_line() +
                gglot2::theme_bw() +
                gglot2::ylim(0, 1) +
                gglot2::labs(x="Sample", y="Proportion in cluster") +
                gglot2::scale_colour_discrete("Cluster", guide=F, drop=F)
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
