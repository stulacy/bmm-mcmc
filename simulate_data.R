# Firstly create a dataset with no noise, i.e. all variables have a different
# theta for each cluster
set.seed(17)
N <- 1000
P <- 5

theta_actual <- matrix(c(0.7, 0.8, 0.2, 0.1, 0.1,
                         0.3, 0.5, 0.9, 0.8, 0.6,
                         0.1, 0.2, 0.5, 0.4, 0.9),
                       nrow=3, ncol=5, byrow = T)
cluster_ratios <- c(0.6, 0.2, 0.2)

mat_list <- lapply(seq_along(cluster_ratios), function(i) {
    n <- round(N * cluster_ratios[i])
    sapply(theta_actual[i,], function(p) rbinom(n, 1, p))
})

mat <- do.call('rbind', mat_list)
saveRDS(mat, "data/K3_N1000_P5_clean.rds")


set.seed(17)
N <- 100
P <- 5

theta_actual <- matrix(c(0.7, 0.8, 0.2, 0.1, 0.1,
                         0.2, 0.2, 0.9, 0.8, 0.6),
                       nrow=2, ncol=5, byrow = T)
cluster_ratios <- c(0.7, 0.3)

mat_list <- lapply(seq_along(cluster_ratios), function(i) {
    n <- round(N * cluster_ratios[i])
    sapply(theta_actual[i,], function(p) rbinom(n, 1, p))
})

mat <- do.call('rbind', mat_list)
saveRDS(mat, "data/K2_N100_P5_clean.rds")

set.seed(17)
N <- 1000
P <- 5

theta_actual <- matrix(c(0.7, 0.8, 0.2, 0.1, 0.1,
                         0.2, 0.2, 0.9, 0.8, 0.6),
                       nrow=2, ncol=5, byrow = T)
cluster_ratios <- c(0.7, 0.3)

mat_list <- lapply(seq_along(cluster_ratios), function(i) {
    n <- round(N * cluster_ratios[i])
    sapply(theta_actual[i,], function(p) rbinom(n, 1, p))
})

mat <- do.call('rbind', mat_list)
saveRDS(mat, "data/K2_N1000_P5_clean.rds")