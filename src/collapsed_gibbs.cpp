// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

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
                         int burnin,
                         bool debug) {

    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;

    // Random sample for first row
    arma::Mat<int> z_out(nsamples, N);
    z_out.row(0) = as<arma::Row<int>>(initialK);
    arma::cube thetas(K, P, nsamples);
    arma::cube probs_out(nsamples, N, K);
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
                    LHS = log(Nk + (alpha/K)) - log(N - 1 + alpha);
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

            // Save probabilities so can fix label switching after sampling
            probs_out.tube(j, i) = as<arma::vec>(probs);

            // Sample z_i using R multinom as slightly more efficient
            // this returns a ones hot encoded binary that needs converting to int
            rmultinom(1, probs.begin(), K, this_z.begin());

            if (debug) Rcout << "probs_sum: " << probs_sum << "\n";
            if (debug) Rcout << "Normalised probs: " << probs << "\n";
            if (debug) Rcout << "sampled z: " << this_z << "\n";

            // Save cluster assignment as label in 1...K rather than one hot encoded binary
            for (int k = 0; k < K; ++k) {
                if (this_z(k) == 1) {
                    z_out(j, i) = k + 1;
                    clusters[k].push_back(i);
                    continue;
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
            }
        }
    }
    
    List ret;
    ret["z"] = z_out.rows(burnin, nsamples-1);
    // Due to how subcubes work (i.e. when have sliced them), can't just
    // return directly, but instead need to force as cube
    arma::cube thetas_post = thetas.tail_slices(nsamples - burnin);
    arma::cube probs_post = probs_out.rows(burnin, nsamples-1);
    ret["theta"] = thetas_post;
    ret["probabilities"] = probs_post;
    return ret;
}

