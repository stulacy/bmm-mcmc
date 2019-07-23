// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

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
                         bool debug) {

    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;

    // Random sample for first row
    arma::Mat<int> allocations(nsamples, N);
    allocations.row(0) = as<arma::Row<int>>(initialK);
    arma::cube thetas(K, P, nsamples);
    int curr_cluster;

    // Vector for each cluster to track members
    // Helps to quickly identify Ck
    std::vector< std::vector< int > > clusters(K);
    for (int i=0; i < N; ++i) {
        clusters[initialK(i)].push_back(i);
    }

    std::vector <int> Ck;
    if (debug) Rcout << "Initial K: " << initialK << "\n";
    if (debug) Rcout << "Initial K: " << initialK << "\n";
    if (debug) Rcout << "N: " << N << "\n";


    // At each sample, for each person:
    for (int j=1; j < nsamples; ++j) {
        Rcout << "Sample " << j+1 << "\n";
        for (int i = 0; i < N; ++i) {
            if (debug) Rcout << "Individual " << i << "\n";
            // Drop that person from current cluster
            curr_cluster = allocations(j-1, i);
            clusters[curr_cluster].erase(std::remove(clusters[curr_cluster].begin(),
                                                     clusters[curr_cluster].end(),
                                                     i),
                                         clusters[curr_cluster].end());

            // For k in 1:K calculate probabilities of being in k by use of the same equation as before
            // Firstly identify set of patients in this cluster
            std::vector <double> probs(K);
            double probs_sum=0;
            if (debug) Rcout << "Calculating probabilities\n";
            for (int k=0; k < K; ++k) {
                if (debug) Rcout << "k: " << k << "\n";
                Ck = clusters[k];
                int Nk = Ck.size();
                double LHS = log(Nk) - log(N - 1 + alpha);
                if (debug) Rcout << "Nk: " << Nk << "\n";
                if (debug) Rcout << "LHS: " << LHS << "\n";
                double logLH = 0;
                for (int d=0; d < P; ++d) {
                    if (debug) Rcout << "d: " << d << "\n";
                    int sum_xd = 0;
                    for (int c : Ck) {
                        sum_xd += df_arma(c, d);
                    }

                    int xnd = df_arma(i, d);
                    double left = xnd * log(beta + sum_xd);
                    double right = (1-xnd) * log(gamma + Nk - sum_xd);
                    double denom = log(beta + gamma + Nk);
                    double full = left + right - denom;

                    if (debug) Rcout << "sum_xd: " << sum_xd << "\n";
                    if (debug) Rcout << "left: " << left << "\n";
                    if (debug) Rcout << "right: " << right << "\n";
                    if (debug) Rcout << "denom: " << denom << "\n";
                    if (debug) Rcout << "full: " << full << "\n";

                    logLH += full;
                }
                double dummy = exp(LHS + logLH);
                probs_sum += dummy;
                if (debug) Rcout << "dummy: " << dummy << "\n";
                probs[k] = dummy;
            }

            if (debug) {
                Rcout << "Raw probs: \t";
                for (auto it : probs) {
                    Rcout << it << " ";
                }
                Rcout << "\n ";
            }

            // Normalise and form K+1 cluster labels to sample from
            IntegerVector choices(K);
            NumericVector probs_norm(K);
            for (int k=0; k < K; ++k) {
                choices(k) = k;
                probs_norm(k) = probs[k] / probs_sum;
            }

            if (debug) Rcout << "probs_sum: " << probs_sum << "\n";
            if (debug) Rcout << "Normalised probs: " << probs_norm << "\n";

            // Sample z_i from it (using R's sample function rather than rmultinom)
            int ret = RcppArmadillo::sample(choices, 1, false, probs_norm)(0);
            if (debug) Rcout << "Sample: " << ret << "\n";

            clusters[ret].push_back(i);
            allocations(j, i) = ret;
        }

        // Estimate thetas
        for (int k=0; k < K; ++k) {
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
    }

    List ret;
    ret["z"] = allocations;
    ret["theta"] = thetas;
    return ret;
}

