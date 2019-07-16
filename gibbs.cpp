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

// Using Algorithm 2 from
// https://www4.stat.ncsu.edu/~wilson/prelim/Review1.pdf
// [[Rcpp::export]]
IntegerMatrix gibbs_cpp(IntegerMatrix df,
                        IntegerVector initialPi,
                        NumericMatrix initialTheta,
                        int nsamples,
                        int K,
                        double alpha,
                        double beta,
                        double gamma,
                        bool verbose) {

    int N = df.nrow();
    int P = df.ncol();


    // Random sample for first row
    NumericMatrix pi_sampled(nsamples, K);
    cube z_sampled(N, K, nsamples, int);
    cube theta_sampled(K, P, nsamples, double);

    int N_nk, sum_x, Znk;
    double frac, num_1, num_2, loglh, denom, contrib, cum_probs, dummy;
    pi_sampled(0, _) = initialPi;
    theta_sampled.slice(0) = initialTheta;
    NumericMatrix thisTheta;
    IntegerVector ck(K);
    IntegerMatrix Vkd(K, D);

    for (int j=1; j < nsamples; j++) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";

        thisTheta = theta_sample.slice(j-1);

        for (int i=0; i < N; i++) {

            //Rcpp::Rcout << "i: " << i << "\n";
            std::vector<double> probs(K);
            NumericVector s(K);
            IntegerVector this_z(K);

            cum_probs = 0;
            loglh = 0;

            // Calculate data likelihood
            for (int k=1; k <= K; k++) {
                for (int d=0; d < P; d++) {
                    loglh += df(i, d) * log(theta(k, d)) + (1 - df(i, d)) * log(1 - theta(k, d));
                }

                // Then calculate probabilities per cluster
                dummy = exp(pi_sampled(j-1, k-1) + loglh);
                s[k-1] = dummy;
                cum_probs += dummy;
            }

            // Then normalise
            for (int p = 0; p < K; p++) {
                s[p] /= cum_probs;
            }

            // Now can draw labels
            rmultinom(1, s.begin(), K, this_z.begin());
            z_sampled.slice(j)(i, _) = this_z;
        }

        // Now sample pi and thetas
        // Now calculate number of patients in each cluster and sum of data points as before
        for (int k = 0; k < K; k++) {
            ck[k] = 0;
            for (int i=0; i < N; i++) {
                Vkd(k, d) = 0;
                Znk = z_sampled.slice(j)(i, k);
                ck[k] += Znk;
                for (int d = 0; d < N; d++) {
                    Vkd(k, d) += Znk * df(i, d);

                }



            }

        }








                //Rcpp::Rcout << "k: " << k << "\n";
                // Setup this set and its cardinality
                N_nk = ck.size();
                //Rcpp::Rcout << "N_nk: " << N_nk << "\n";

                // Second half is loglh
                loglh = 0;
                for (int d=0; d < P; d++) {
                    sum_x = calc_sumx(df, ck, d);
                }
            }

            if (verbose) {
                Rcpp::Rcout << "raw probs: ";
                for (double foo : probs) {
                     Rcpp::Rcout << foo << " | ";
                }
                Rcpp::Rcout << "\n";
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

