// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <stdlib.h>
#include <queue>
#include "stephens.h"
#include "utils.h"

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
                            double a,
                            double b,
                            int burnin,
                            bool relabel,
                            int burnrelabel,
                            int maxK,
                            bool debug) {

    // Setup int K giving current number of clusters, initialised to number of observations
    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;
    int K = 0;
    int curr_cluster;

    if (beta != gamma) {
        Rcpp::stop("Error: sampler currently not implemented for non-symmetric priors on beta and gamma\n");
    }

    // Also want to save allocations along sampler
    arma::Mat<int> allocations(nsamples, N);
    arma::Mat<int> allocations_relabelled(nsamples, N);
    
    // Create 2D vector of cluster membership. Vector is K length with each entry being
    // 1D vector detailing the observations residing in that cluster
    std::vector< std::vector< int > > clusters(maxK);
    // a 1D vector of ints max size N, detailing the clusters that are in use
    // i.e. which indices of clusters are not empty
    std::vector<int> used_clusters;
    std::priority_queue<int, std::vector<int>, std::greater<int> > unused_clusters; 
    for (int i = 0; i < N; ++i) {
        unused_clusters.push(i);
    }

    if (debug) print_clusters(clusters);

    // Can calculate probability of new cluster outside of loop as is constant
    // due to use of symmetric priors on theta with Beta(beta, gamma) beta == gamma
    double RHS_newk = P * (log(beta) - log(beta + gamma));
    
    int xnd, sum_xd, Nk, smallest_cluster, smallest_size, this_cluster_size;
    double left, right, LHS, denom, logLH, max_prob, probs_newk;

    // And want to save probabilities and thetas
    arma::cube thetas(maxK, P, nsamples, arma::fill::zeros);
    arma::cube thetas_relab(maxK, P, nsamples, arma::fill::zeros);

    std::vector <int> Ck;
    double b_eps, pi, pi1, pi2;
    arma::vec alpha_sampled(nsamples);
    alpha_sampled(0) = alpha;
    double alpha_new, foobar, sumprob, left_denom;
    arma::cube probs_out(N, maxK, burnrelabel, arma::fill::zeros);
    arma::mat probs_sample(N, maxK, arma::fill::zeros);
    arma::mat Q;
    std::pair<arma::Row<int>, arma::mat> stephens_out;
    arma::Row<int> permutations_sample(maxK);
    arma::Mat<int> permutations(nsamples - burnin, maxK);
        
    for (int j=1; j < nsamples; ++j) {
        Rcout << "Sample " << j+1 << "\tK: " << K << "\n";

        // Constant value on denominator of left hand fraction updates with alpha each sample
        left_denom = log(N - 1 + alpha_sampled(j-1));

        // Probability of creating a new cluster using Eq 25 from Maartens only updates with each sample
        // as we're using symmetrical prior on theta
        probs_newk = log(alpha_sampled(j-1)) - left_denom + RHS_newk;

        for (int i = 0; i < N; ++i) {
            
            if (debug) Rcout << "Individual: " << i << "\n";
            
            if (j > 1) {
                // Drop person from current cluster
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
                    unused_clusters.push(curr_cluster);
                    K--;
                }
    
                if (debug) Rcout << "Dropped individual. K: " << K << "\tLength of used_clusters: " << used_clusters.size() << "\n";
            }
    
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
                if (Nk == 0) Rcout << "Something's gone wrong... Nk of 0!\n";
                LHS = log(Nk) - left_denom;
                logLH = 0;
                denom = log(beta + gamma + Nk);
                for (int d=0; d < P; ++d) {
                    //if (debug) Rcout << "d: " << d << "\n";
                    sum_xd = 0;
                    for (int c : Ck) {
                        sum_xd += df_arma(c, d);
                    }

                    xnd = df_arma(i, d);
                    left = xnd * log(beta + sum_xd);
                    right = (1-xnd) * log(gamma + Nk - sum_xd);
                    logLH += left + right - denom;
                }
                probs(k) = LHS + logLH;
                choices(k) = used_clusters[k];
            }

            // Add on probability and label for new K, which will label as an unused
            // cluster in the N dimensions, and hence is available
            if (unused_clusters.size() == 0) {
                Rcpp::stop("Error: have no free clusters, need to create one.");
            }
            int new_cluster = unused_clusters.top();
            choices(K) = new_cluster;
            probs(K) = probs_newk;

            // Calculate exponentiated probs using exponentiate-normalise trick
            max_prob = max(probs);
            if (debug) Rcout << "Max prob: " << max_prob << "\n";
            sumprob=0;
            for (int k = 0; k <= K; ++k) {
                foobar = exp(probs(k) - max_prob);
                probs_norm(k) = foobar;
                sumprob += foobar;
            }

            // Finish softmax exp-normalise 
            for (int k=0; k <= K; ++k) {
                probs_norm(k) /= sumprob;
            }
            
            // Save initial probabilities so can do batch Stephens
            // to generate initial Q
            if (relabel) {
                if (j < burnin && j >= (burnin - burnrelabel)) {
                    for (int k=0; k <= K; ++k) {
                        probs_out(i, choices(k), j - burnin + burnrelabel) = probs_norm(k);
                    }
                } else if (j >= burnin) {
                    for (int k=0; k <= K; ++k) {
                        probs_sample(i, choices(k)) = probs_norm(k);
                    }
                }
            }
            
            if (debug) Rcout << "Raw probs (" << probs.size() << ") :" << probs << "\n";
            if (debug) Rcout << "Normalised probs (" << probs_norm.size() << ") :" << probs_norm << "\n";
            if (debug) Rcout << "Sampling. Length choices: " << choices.size() << "\tLength used_clusters: " << used_clusters.size() << "\tLength unused clusters: " << unused_clusters.size() << "\tK: " << K << "\n";

            // Sample z_i from it (using R's sample function rather than rmultinom)
            int ret = RcppArmadillo::sample(choices, 1, false, probs_norm)(0);
            if (debug) Rcout << "Sampled k: " << ret << "\n";

            // If z_i = K+1 then update list of used and unused clusters
            // But truncate if have enough clusters and add to smallest cluster
            if (ret == new_cluster) {
                if (K < (maxK - 1)) {
                    if (debug) Rcout << "Have created new cluster\n";
                    unused_clusters.pop();
                    used_clusters.push_back(new_cluster);
                    K++;
                } else {
                    if (debug) Rcout << "Already at max K limit so will not create cluster\n";
                    smallest_cluster=0;
                    smallest_size = N+1;
                    for (int k = 0; k < K; ++k) {
                        this_cluster_size = clusters[used_clusters[k]].size();
                        if (this_cluster_size < smallest_size) {
                            smallest_size = this_cluster_size;
                            smallest_cluster = k;
                        }
                        ret = smallest_cluster;
                    }
                }
            }
            clusters[ret].push_back(i);
            allocations(j, i) = ret + 1;
            if (debug) Rcout << "After sampling. Length choices: " << choices.size() << "\tLength used_clusters: " << used_clusters.size() << "\tLength unused clusters: " << unused_clusters.size() << "\tK: " << K << "\n";

            alpha_new = update_alpha(alpha_sampled(j-1), a, b, N, K);
            if (debug) Rcout << "b-log(epsilon): " << b_eps << "\tpi: " << pi << "\tALPHA: " << alpha_new << "\n";
            alpha_sampled(j) = alpha_new;
            
        }  // End for 1:N loop
        
        // To relabel clusters use Stephen's 2000b online algorithm.
        // Firstly need to initialise Q with values taken from a batch formulation
        // over $burnrelabel samples
        if (relabel) {
            if (j == (burnin - 1)) {
                Rcout << "Running Stephens Batch relabelling to identify initial Q values\n";
                Q = my_stephens_batch(probs_out, false);
            } else if (j >= burnin) {
                stephens_out = my_stephens_online(Q, probs_sample, j, false);
                Q = stephens_out.second;
                permutations_sample = stephens_out.first;
                // Relabel Z 
                permutations.row(j-burnin) = permutations_sample;
                for (int i = 0; i < N; ++i) {
                    allocations_relabelled(j, i) = permutations_sample(allocations(j, i) - 1) + 1;
                }
            }
        }
        
        // Estimate thetas
        if (debug) Rcout << "Estimating thetas\n";
        int cluster_index;
        for (int k : used_clusters) {
            Ck = clusters[k];
            int Nk = Ck.size();
            
            for (int d=0; d < P; ++d) {
                int dsum = 0;
                for (int c : Ck) {
                    dsum += df_arma(c, d);
                }
                thetas(k, d, j) = dsum / (double)Nk;
                
                if (relabel && j >= burnin) {
                    thetas_relab(permutations_sample(k), d, j) = dsum / (double)Nk;
                }
            }
        }
        
    }  // End for 1:Nsamples loop
    
    List out;
    arma::cube thetas_post = thetas.tail_slices(nsamples - burnin);
    arma::cube thetas_relabelled = thetas_relab.tail_slices(nsamples - burnin);
    if (relabel) {
        out["z_relabelled"] = allocations_relabelled.tail_rows(nsamples-burnin);
        out["theta_relabelled"] = thetas_relabelled;
    }
    out["z"] = allocations.tail_rows(nsamples-burnin);
    out["theta"] = thetas_post;
    out["alpha"] = alpha_sampled.tail(nsamples-burnin);
    out["permutations"] = permutations;
    return out;
}

