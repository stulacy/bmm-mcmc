// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int calc_sumx(IntegerMatrix data, std::vector<int> ck, int d) {
    int sum_x = 0;
    for (int i : ck) {
        sum_x += data(i, d);
    }
    return sum_x;
}

std::vector<int> get_ck(IntegerMatrix labels, int N, int step, int obs, int k) {
    std::vector<int> out;
    for (int i = 0; i < N; i++) {
        if (i == obs) {
            continue;
        }
        if (labels(step-1, i) == k) {
            out.push_back(i);
        }
    }
    return out;
}

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
IntegerMatrix collapsed_gibbs_cpp(IntegerMatrix df,
                                  IntegerVector initialK,
                                  int nsamples,
                                  int K,
                                  double alpha,
                                  double beta,
                                  double gamma,
                                  bool verbose) {

    int N = df.nrow();
    int P = df.ncol();

    // Random sample for first row
    IntegerMatrix allocations(nsamples, N);
    int N_nk, sum_x;
    std::vector<int> ck;
    double frac, num_1, num_2, loglh, denom, contrib, cum_probs, dummy;
    allocations(0, _) = initialK;
    for (int j=1; j < nsamples; j++) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";
        // Set params
        for (int i=0; i < N; i++) {
        //Rcpp::Rcout << "i: " << i << "\n";
            std::vector<double> probs(K);
            cum_probs = 0;
            for (int k=1; k <= K; k++) {
                //Rcpp::Rcout << "k: " << k << "\n";
                // Setup this set and its cardinality
                ck = get_ck(allocations, N, j, i, k);
                N_nk = ck.size();
                //Rcpp::Rcout << "N_nk: " << N_nk << "\n";

                // First half
                frac = log(N_nk + alpha/K) - log(N - 1 + alpha);
                if (verbose) Rcpp::Rcout << "Frac: " << frac << "\n";

                // Second half is loglh
                loglh = 0;
                for (int d=0; d < P; d++) {
                    sum_x = calc_sumx(df, ck, d);
                    //Rcpp::Rcout << "sum_x: " << sum_x << "\n";

                    num_1 = df(i,d) * log(beta + sum_x);
                    //Rcpp::Rcout << "num_1: " << num_1 << "\n";
                    num_2 = (1-df(i,d)) * log(gamma + N_nk - sum_x);
                    //Rcpp::Rcout << "num_2: " << num_2 << "\n";
                    denom = log(beta + gamma + N_nk);
                    //Rcpp::Rcout << "denom: " << denom << "\n";
                    contrib = num_1 + num_2 - denom;
                    //Rcpp::Rcout << "contrib: " << contrib << "\n";
                    loglh = loglh + contrib;
                }
                if (verbose) Rcpp::Rcout << "loglh: " << loglh << "\n";
                dummy = exp(frac+loglh);
                probs[k-1] = dummy;
                cum_probs += dummy;
            }
            if (verbose) {
                Rcpp::Rcout << "raw probs: ";
                for (double foo : probs) {
                     Rcpp::Rcout << foo << " | ";
                }
                Rcpp::Rcout << "\n";
            }
            for (int p = 0; p < K; p++) {
                probs[p] /= cum_probs;
            }
            if (verbose) {
                Rcpp::Rcout << "normalised probs: ";
                for (double foo : probs) {
                     Rcpp::Rcout << foo << " | ";
                }
                Rcpp::Rcout << "\n";
            }

            // Can't seem to pass in seq_len to first arg
            // so need to manually form IntegerVector
            // to sample form
            IntegerVector vals(K);
            for (int i = 0; i < K; i++) {
                vals(i) = i + 1;
            }
            allocations(j, i) = Rcpp::RcppArmadillo::sample(vals, 1, false, wrap(probs))(0);
        }
    }
    return allocations;
}

