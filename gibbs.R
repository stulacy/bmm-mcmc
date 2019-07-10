# Runs Gibbs sampling
# Starting with fixed K
# Then will buildup to random K using Dirichlet prior
# Then will finalise with infinite mixture using Dirichlet Process prior

df <- readRDS("data/N1000_P5_clean.rds")

N <- nrow(df)
P <- ncol(df)

# Since we simulated the data we know what K actually is
K_actual <- 3


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

# So HOW DO WE SAMPLE FROM THIS?
# Set initial values

nsamples <- 1000
K <- 3

alpha <- 1
beta <- 1
gamma <- 1

# Random sample for first row
allocations <- matrix(0, nrow=nsamples, ncol=N)
allocations[1, ] <- sample(1:K, N, replace=T)
for (j in 2:nsamples) {
    cat(paste0("Sample: ", j, "\n"))
    # Set params
    for (i in 1:N) {
        probs <- rep(0, K)
        for (k in 1:K) {
            # All points in current cluster except current data point. HOW DOES THIS DIFFER FROM N_nk?
            ck <- allocations[j-1, -i] == k
            N_nk <- sum(ck)

            frac <- (N_nk * alpha/K) / (N - 1 + alpha)

            loglh <- 0
            for (d in 1:P) {
                sum_x <- sum(df[-i, ][ck, d])

                num_1 <- (beta + sum_x)**(df[i, d])
                num_2 <- (gamma + N_nk + sum_x)**(1 - df[i, d])
                contrib <- log(num_1) + log(num_2) - log(beta + gamma + N_nk)
                loglh <- loglh + contrib
            }
            probs[k] <- frac * exp(loglh)
            probs <- (probs - min(probs)) / (max(probs) - min(probs))
            #cat(paste0("Probs: ", probs, "\n"))
        }

        allocations[j, i] <- sample(1:k, 1, prob=probs)
    }
}
