// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"
#include "stephens.h"

using namespace Rcpp;

// Collapsed Gibbs sampler
// For each i in 1:N
//  Sample z_i given (theta, z_(-i), x)
// For k in 1:K
//  Sample theta_k given (theta_(-k), z, x)
//
// Tu in Section 3.23 says needs to sample p(z | ...), p(pi | ...), and p(theta | ...)
// But if use collapsed Gibbs then can sample p(z) directly
// https://people.eecs.berkeley.edu/~stephentu/writeups/mixturemodels.pdf
//
// Eq 21 from Van Maarten also suggests can just use collapsed Gibbs sampler over z
// this equation looks to be the same as Tu's version.
// https://pdfs.semanticscholar.org/525c/ff658d34ae4d47e84b8ec4ede3ce6c561afc.pdf
//
// And this looks to be the same as 3.5 from Neal
// http://www.stat.columbia.edu/npbayes/papers/neal_sampling.pdf
// [[Rcpp::export]]
List collapsed_gibbs_cpp(IntegerMatrix df,
                         IntegerVector initialK,
                         int nsamples,
                         int K,
                         double alpha,
                         double beta,
                         double gamma,
                         double a, 
                         double b,
                         int burnin,
                         bool relabel,
                         int burnrelabel,
                         bool debug) {

    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;

    // Random sample for first row
    arma::Mat<int> z_out(nsamples, N);
    arma::Mat<int> z_out_relabelled(nsamples, N);
    z_out.row(0) = as<arma::Row<int>>(initialK);
    arma::cube thetas(K, P, nsamples);
    arma::cube thetas_relab(K, P, nsamples);
    arma::vec alpha_sampled(nsamples);
    if (alpha == 0) {
        alpha_sampled(0) = 1;
    } else {
        alpha_sampled.fill(alpha);
    }

    int curr_cluster;

    // Vector for each cluster to track members
    // Helps to quickly identify Ck
    std::vector< std::vector< int > > clusters(K);
    for (int i=0; i < N; ++i) {
        clusters[initialK(i)-1].push_back(i);
    }

    std::vector <int> Ck;
    if (debug) Rcout << "Initial K: " << initialK << "\n";
    if (debug) Rcout << "N: " << N << "\n";

    NumericVector probs(K);
    IntegerVector this_z(K);
    double probs_sum;
    int Nk, sum_xd, xnd, dsum;
    double LHS, logLH, left, right, denom, full, dummy;
    
    // Data structures for relabelling
    arma::cube probs_out(N, K, burnrelabel, arma::fill::zeros);
    arma::mat probs_sample(N, K, arma::fill::zeros);
    arma::mat Q;
    std::pair<arma::Row<int>, arma::mat> stephens_out;
    arma::Row<int> permutations_sample(K);
    arma::Mat<int> permutations(nsamples - burnin, K);

    // At each sample, for each person:
    for (int j=1; j < nsamples; ++j) {
        Rcout << "Sample " << j+1 << "\n";
        for (int i = 0; i < N; ++i) {
            if (debug) Rcout << "Individual " << i << "\n";
            // Drop that person from current cluster
            curr_cluster = z_out(j-1, i) - 1;
            clusters[curr_cluster].erase(std::remove(clusters[curr_cluster].begin(),
                                                     clusters[curr_cluster].end(),
                                                     i),
                                         clusters[curr_cluster].end());

            // For k in 1:K calculate probabilities of being in k by use of the same equation as before
            // Firstly identify set of patients in this cluster
            probs_sum=0;
            if (debug) Rcout << "Calculating probabilities\n";
            for (int k=0; k < K; ++k) {
                if (debug) Rcout << "k: " << k << "\n";
                Ck = clusters[k];
                Nk = Ck.size();
                
                if (Nk > 0) {
                    LHS = log(Nk + (alpha_sampled(j-1)/K)) - log(N - 1 + alpha_sampled(j-1));
                    if (debug) Rcout << "Nk: " << Nk << "\n";
                    if (debug) Rcout << "LHS: " << LHS << "\n";
                    logLH = 0;
                    for (int d=0; d < P; ++d) {
                        if (debug) Rcout << "d: " << d << "\n";
                        sum_xd = 0;
                        for (int c : Ck) {
                            sum_xd += df_arma(c, d);
                        }
    
                        xnd = df_arma(i, d);
                        left = xnd * log(beta + sum_xd);
                        right = (1-xnd) * log(gamma + Nk - sum_xd);
                        denom = log(beta + gamma + Nk);
                        full = left + right - denom;
    
                        if (debug) Rcout << "sum_xd: " << sum_xd << "\n";
                        if (debug) Rcout << "left: " << left << "\n";
                        if (debug) Rcout << "right: " << right << "\n";
                        if (debug) Rcout << "denom: " << denom << "\n";
                        if (debug) Rcout << "full: " << full << "\n";
    
                        logLH += full;
                    }
                    dummy = exp(LHS + logLH);
                } else {
                    dummy = 0;
                }
                probs_sum += dummy;
                if (debug) Rcout << "dummy: " << dummy << "\n";
                probs(k) = dummy;
            }

            if (debug) {
                Rcout << "Raw probs: \t";
                for (auto it : probs) {
                    Rcout << it << " ";
                }
                Rcout << "\n ";
            }

            // Normalise probs
            for (int k=0; k < K; ++k) {
                probs(k) /= probs_sum;
            }

            // Sample z_i using R multinom as slightly more efficient
            // this returns a ones hot encoded binary that needs converting to int
            rmultinom(1, probs.begin(), K, this_z.begin());

            if (debug) Rcout << "probs_sum: " << probs_sum << "\n";
            if (debug) Rcout << "Normalised probs: " << probs << "\n";
            if (debug) Rcout << "sampled z: " << this_z << "\n";
            
            // Save initial probabilities so can do batch Stephens
            // to generate initial Q
            if (relabel) {
                if (j < burnin && j >= (burnin - burnrelabel)) {
                    for (int k=0; k < K; ++k) {
                        probs_out(i, k, j - burnin + burnrelabel) = probs(k);
                    }
                } else if (j >= burnin) {
                    for (int k=0; k < K; ++k) {
                        probs_sample(i, k) = probs(k);
                    }
                }
            }

            // Save cluster assignment as label in 1...K rather than one hot encoded binary
            for (int k = 0; k < K; ++k) {
                if (this_z(k) == 1) {
                    z_out(j, i) = k + 1;
                    clusters[k].push_back(i);
                    continue;
                }
            }
        }  // End 1:N loop
        
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
                    z_out_relabelled(j, i) = permutations_sample(z_out(j, i) - 1) + 1;
                }
            }
        }
        

        // Estimate thetas
        for (int k=0; k < K; ++k) {
            Ck = clusters[k];
            Nk = Ck.size();
            for (int d=0; d < P; ++d) {
                dsum = 0;
                for (int c : Ck) {
                    dsum += df_arma(c, d);
                }
                if (debug) Rcout << "Theta k: " << k << "\td: " << d << "\tdsum: " << dsum << "\tNk: " << Nk << "\tdsum / NK: " << dsum / (double)Nk << "\n";
                thetas(k, d, j) = dsum / (double)Nk;
                if (relabel && j >= burnin) {
                    thetas_relab(permutations_sample(k), d, j) = dsum / (double)Nk;
                }
            }
        }
        
        // Update alpha
        if (alpha == 0) {
            alpha_sampled(j) = update_alpha(alpha_sampled(j-1), a, b, N, K);
        }
    }
    
    // Due to how subcubes work (i.e. when have sliced them), can't just
    // return directly, but instead need to force as cube
    List ret;
    arma::cube thetas_post = thetas.tail_slices(nsamples - burnin);
    ret["alpha"] = alpha_sampled.tail_rows(nsamples-burnin);
    ret["permutations"] = permutations;
    if (relabel) {
        arma::cube thetas_relabelled = thetas_relab.tail_slices(nsamples - burnin);
        ret["z"] = z_out_relabelled.tail_rows(nsamples-burnin);
        ret["theta"] = thetas_relabelled;
        ret["z_original"] = z_out.tail_rows(nsamples-burnin);
        ret["theta_original"] = thetas_post;
    } else {
        ret["z"] = z_out.tail_rows(nsamples-burnin);
        ret["theta"] = thetas_post;
    }
    return ret;
}

