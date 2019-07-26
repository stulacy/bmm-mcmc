// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <stdlib.h>

using namespace Rcpp;

void print_clusters(std::vector < std::vector < int > > clusters) {
    for (unsigned int i=0; i < clusters.size(); i++) {
        Rcout << "Cluster: " << i << "\n";
        for (auto it : clusters[i]) {
            Rcout << "Individual: " << it << "\n";
        }
    }
}


// Collapsed Gibbs sampler using Dirichlet Process Prior on cluster weights
// Using Algorithm 3 from Neal as basis
// http://www.stat.columbia.edu/npbayes/papers/neal_sampling.pdf
// Where the two integrals have been expanded by Maartens in eq 23 below
// https://pdfs.semanticscholar.org/525c/ff658d34ae4d47e84b8ec4ede3ce6c561afc.pdf
// [[Rcpp::export]]
List collapsed_gibbs_dp_cpp(IntegerMatrix df,
                            int nsamples,
                            double alpha,
                            double beta,
                            double gamma,
                            bool debug) {

    // Setup int K giving current number of clusters, initialised to number of observations
    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;
    int K = N;
    int curr_cluster;

    if (beta != gamma) {
        Rcpp::stop("Error: sampler currently not implemented for non-symmetric priors on beta and gamma\n");
    }

    // Also want to save allocations along sampler
    arma::Mat<int> allocations(nsamples, N);

    // Create 2D vector of cluster membership. Vector is K length with each entry being
    // 1D vector detailing the observations residing in that cluster
    std::vector< std::vector< int > > clusters(N);
    // a 1D vector of ints max size N, detailing the clusters that are in use
    // i.e. which indices of clusters are not empty
    std::vector<int> used_clusters;
    std::vector<int> unused_clusters;
    for (int i = 0; i < N; ++i) {
        clusters[i] = std::vector<int> {i};
        allocations(0, i) = i + 1;
        used_clusters.push_back(i);
    }

    if (debug) print_clusters(clusters);

    // Can calculate probability of new cluster outside of loop as is constant
    // due to use of symmetric priors on theta with Beta(beta, gamma) beta == gamma
    double RHS_newk = P * (log(beta) - log(beta + gamma));
    double LHS_newk;
    double dummy_newk;


    int xnd, sum_xd, Nk;
    double left, right, dummy, LHS, logLH, max_prob;

    arma::cube thetas(K, P, nsamples, arma::fill::zeros);
    std::vector <int> Ck;
    double a=1;
    double b=1;
    double epsilon, pi1, pi2, pi;
    NumericVector alpha_sampled(nsamples);
    alpha_sampled(0) = alpha;
    double alpha_new, foobar, sumprob;

    // At each sample, for each person:
    for (int j=1; j < nsamples; ++j) {
        Rcout << "Sample " << j+1 << "\tK: " << K << "\n";

        for (int i = 0; i < N; ++i) {

            // Drop that person from current cluster
            curr_cluster = allocations(j-1, i) - 1;
            clusters[curr_cluster].erase(std::remove(clusters[curr_cluster].begin(),
                                                     clusters[curr_cluster].end(),
                                                     i),
                                         clusters[curr_cluster].end());

            // If this makes it empty, then remove it from list of used clusters and decrement K
            if (clusters[curr_cluster].size() == 0) {

                used_clusters.erase(std::remove(used_clusters.begin(),
                                                used_clusters.end(),
                                                curr_cluster),
                                    used_clusters.end());
                unused_clusters.push_back(curr_cluster);
                K--;
            }

            if (debug) Rcout << "Dropped individual. K: " << K << "\tLength of used_clusters: " << used_clusters.size() << "\n";

            //print_clusters(clusters);

            // For k in 1:K calculate probabilities of being in k by use of the same equation as before
            // Firstly identify set of patients in this cluster
            NumericVector probs(K+1);
            IntegerVector choices(K+1);
            NumericVector probs_norm(K+1);

            //if (debug) Rcout << "Calculating probabilities for known K\n";
            for (int k = 0; k < K; ++k) {
                //if (debug) Rcout << "k: " << k << "\n";
                Ck = clusters[used_clusters[k]];
                Nk = Ck.size();
                LHS = log(Nk) - log(N - 1 + alpha_sampled(j-1));
                //if (debug) Rcout << "Nk: " << Nk << "\n";
                //if (debug) Rcout << "LHS: " << LHS << "\n";
                logLH = 0;
                for (int d=0; d < P; ++d) {
                    //if (debug) Rcout << "d: " << d << "\n";
                    sum_xd = 0;
                    for (int c : Ck) {
                        sum_xd += df_arma(c, d);
                    }

                    xnd = df_arma(i, d);
                    left = xnd * log(beta + sum_xd);
                    right = (1-xnd) * log(gamma + Nk - sum_xd);
                    logLH += left + right - log(beta + gamma + Nk);
                }
                dummy = LHS + logLH;
                if (debug) Rcout << "dummy: " << dummy << "\n";
                probs(k) = dummy;
            }

            // Add on probability of creating a new cluster using Eq 25 from Maartens
            LHS_newk = log(alpha_sampled(j-1)) - log(N - 1 + alpha_sampled(j-1));
            dummy_newk = LHS_newk + RHS_newk;
            probs(K) = dummy_newk;

            // Calculate exponentiated probs using exponentiate-normalise trick
            max_prob = max(probs);
            if (debug) Rcout << "Max prob: " << max_prob << "\n";
            sumprob=0;
            for (int k = 0; k <= K; ++k) {
                foobar = exp(probs(k) - max_prob);
                probs_norm(k) = foobar;
                sumprob += foobar;
            }

            // Finish softmax exp-normalise and form K cluster labels to sample from
            for (int k=0; k < K; ++k) {
                choices(k) = used_clusters[k];
                probs_norm(k) /= sumprob;
            }
            probs_norm(K) /= sumprob;

            // Add on probability and label for new K, which will label as an unused
            // cluster in the N dimensions, and hence is available
            if (unused_clusters.size() == 0) {
                Rcpp::stop("Error: have no free clusters, need to create one.");
            }
            int new_cluster = unused_clusters.back();
            choices(K) = new_cluster;

            //if (debug) Rcout << "probs_sum: " << probs_sum << "\n";
            if (debug) Rcout << "Raw probs (" << probs.size() << ") :" << probs << "\n";
            if (debug) Rcout << "Normalised probs (" << probs_norm.size() << ") :" << probs_norm << "\n";
            if (debug) Rcout << "Sampling. Length choices: " << choices.size() << "\tLength used_clusters: " << used_clusters.size() << "\tLength unused clusters: " << unused_clusters.size() << "\tK: " << K << "\n";

            // Sample z_i from it (using R's sample function rather than rmultinom)
            int ret = RcppArmadillo::sample(choices, 1, false, probs_norm)(0);
            if (debug) Rcout << "Sample: " << ret << "\n";

            // If z_i = K+1 then update list of used and unused clusters
            if (ret == new_cluster) {
                if (debug) Rcout << "Have created new cluster\n";
                unused_clusters.pop_back();
                used_clusters.push_back(new_cluster);
                K++;
            }
            clusters[ret].push_back(i);
            allocations(j, i) = ret + 1;
            if (debug) Rcout << "After sampling. Length choices: " << choices.size() << "\tLength used_clusters: " << used_clusters.size() << "\tLength unused clusters: " << unused_clusters.size() << "\tK: " << K << "\n";
        }

        // Estimate thetas
        for (int k : used_clusters) {
            Ck = clusters[k];
            int Nk = Ck.size();
            for (int d=0; d < P; ++d) {
                int dsum = 0;
                for (int c : Ck) {
                    dsum += df_arma(c, d);
                }
                if (debug) Rcout << "Theta k: " << k << "\td: " << d << "\tdsum: " << dsum << "\tNk: " << Nk << "\tdsum / NK: " << dsum / (double)Nk << "\n";
                thetas(k, d, j) = dsum / (double)Nk;
            }
        }

        // Sample alpha from Gamma(a, b) using method described by Escobar and West
        // in Section 6
        // https://pdfs.semanticscholar.org/df25/adb36860c1ad9edaac04b8855a2f19e79c5b.pdf
        epsilon = log(R::rbeta(alpha_sampled(j-1)+1, N));
        pi1 = a + K + 1;
        pi2 = N*(b - epsilon);
        pi = pi1 / (pi1 + pi2);
        alpha_new = pi * R::rgamma(a+K, 1/(b-epsilon)) + (1-pi) * R::rgamma(a+K-1, 1/(b-epsilon));
        if (debug) Rcout << "log(epsilon): " << epsilon << "\tpi1: " << pi1 << "\tpi2: " << pi2 << "\tpi: " << pi << "\tALPHA: " << alpha_new << "\n";
        alpha_sampled(j) = alpha_new;

    }
    List out;
    out["z"] = allocations;
    out["theta"] = thetas;
    out["alpha"] = alpha_sampled;
    return out;
}

