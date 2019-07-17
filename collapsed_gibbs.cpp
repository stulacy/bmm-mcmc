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

int calc_sumx(arma::Mat<int> data, std::vector<int> ck, int d, int curr_sample) {
    int sum_x = 0;
    for (int i : ck) {
        if (i == curr_sample) continue;
        sum_x += data(i, d);
    }
    return sum_x;
}

std::vector<int> get_ck(IntegerMatrix labels, int N, int step, int obs, int k) {
    std::vector<int> out;
    for (int i = 0; i < N; ++i) {
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
List collapsed_gibbs_cpp(IntegerMatrix df,
                         IntegerMatrix initialK,
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
    arma::Cube<int> allocations(N, K, nsamples);
    allocations.slice(0) = as<arma::Mat<int>>(initialK);
    arma::cube thetas(K, P, nsamples);
    
    int N_nk, Nk, sum_x, dpoint;
    std::vector<int> this_ck;
    double frac, num_1, num_2, loglh, denom, contrib, cum_probs, dummy;
    for (int j=1; j < nsamples; ++j) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";
        // Set params
        for (int i=0; i < N; ++i) {
            NumericVector probs(K);
            cum_probs = 0;
            arma::uvec ck_inds = arma::linspace<arma::uvec>(0, N-1, N);
            ck_inds.shed_row(i);
            Rcout << "ck inds size: " << arma::size(ck_inds) << "\n";
            
            for (int k=0; k < K; ++k) {
                arma::Col<int> Ck = allocations.slice(j-1).rows(ck_inds).col(k);
                Rcout << "Ck size: " << arma::size(Ck) << "\n";
                Nk = arma::sum(Ck);
                N_nk = Nk - allocations(i, k, j-1);
                if (debug) Rcout << "Nk: " << Nk << "\tN_nk: " << N_nk << "\n";

                // First half
                frac = log(N_nk + alpha/K) - log(N - 1 + alpha);
                if (debug) Rcpp::Rcout << "Frac: " << frac << "\n";

                // Second half is loglh
                loglh = 0;
                
                // TODO Correctly remove current data point from Ck
                // TODO Vectorise this
                for (int d=0; d < P; ++d) {
                    // Ck = (N,)
                    // data[, d] = (N)
                    
                    
                    
                    arma::Col<int> out = df_arma.col(d).t() * Ck;
                    sum_x = out(0,0);
                    dpoint = df_arma(i,d);
                    if (debug) Rcout << "d: " << d << "\tsum_x: " << sum_x << "\tNk: " << Nk << "\n";
                    thetas(k, d, j) = sum_x / Nk;

                    num_1 = dpoint * log(beta + sum_x);
                    num_2 = (1-dpoint) * log(gamma + N_nk - sum_x);
                    denom = log(beta + gamma + N_nk);
                    contrib = num_1 + num_2 - denom;
                    
                    if (debug) Rcout << "num_1: " << num_1;
                    if (debug) Rcout << "\tnum_2: " << num_2;
                    if (debug) Rcout << "\tdenom: " << denom;
                    if (debug) Rcout << "\tcontrib: " << contrib;
                    if (debug) Rcout << "\tcurr loglh: " << loglh << "\n";
                    
                    
                    loglh += contrib;
                }
                if (debug) Rcpp::Rcout << "loglh: " << loglh << "\n";
                dummy = exp(frac+loglh);
                probs(k) = dummy;
                cum_probs += dummy;
            }
            if (debug) {
                Rcpp::Rcout << "raw probs: ";
                for (int k=0; k < K; ++k) {
                     Rcpp::Rcout << probs(k) << " | ";
                }
                Rcpp::Rcout << "\n";
            }
            for (int p = 0; p < K; ++p) {
                probs(p) /= cum_probs;
            }
            if (debug) {
                Rcpp::Rcout << "normalised probs: ";
                for (int k=0; k < K; ++k) {
                     Rcpp::Rcout << probs(k) << " | ";
                }
                Rcpp::Rcout << "\n";
            }

            IntegerVector this_z(K);
            R::rmultinom(1, probs.begin(), K, this_z.begin());
            if (debug) Rcout << "this_z: " << this_z << "\n";
            // TODO Change Col to Row here and can remove transpose?
            allocations.slice(j).row(i) = as<arma::Col<int>>(this_z).t();
            //for (int k=0; k < K; ++k) {
            //    if (this_z(k) == 1) {
            //        if (debug) Rcout << "Adding " << i << " to Ck for k = " << k << "\n";
            //        Ck[j][k].push_back(i);
            //    }
            //}
        }
    }
    
    // Now calculate thetas
    List ret;
    ret["z"] = allocations;
    ret["theta"] = thetas;
    return ret;
}

