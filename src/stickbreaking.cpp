// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Using Algorithm 5.2 from
// http://people.ee.duke.edu/~lcarin/Yuting3.3.06.pdf
// [[Rcpp::export]]
List gibbs_stickbreaking_cpp(IntegerMatrix df,
                             NumericVector initialPi,
                             NumericMatrix initialTheta,
                             int nsamples,
                             int maxK,
                             double alpha,
                             double beta,
                             double gamma,
                             int burnin,
                             bool debug) {

    int N = df.nrow();
    int P = df.ncol();


    // Random sample for first row
    NumericMatrix pi_sampled(nsamples, maxK);
    arma::cube z_sampled(N, maxK, nsamples);
    arma::cube theta_sampled(maxK, P, nsamples);
    arma::Mat<int> z_out(nsamples, N);

    int Znk;
    double loglh, cum_probs, dummy;
    pi_sampled(0, _) = initialPi;
    theta_sampled.slice(0) = as<arma::mat>(initialTheta);
    NumericMatrix thisTheta(maxK,P);
    IntegerVector this_z(maxK);
    NumericVector s(maxK);
    arma::vec dirich_params = arma::zeros(maxK);
    NumericVector this_pi(maxK);
    NumericMatrix theta_row(maxK, P);

    for (int j=1; j < nsamples; ++j) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";

        thisTheta = wrap(theta_sampled.slice(j-1));
        for (int i=0; i < N; ++i) {
            if (debug) Rcout << "Iterating through N to sample z\n";
            cum_probs = 0;

            // Calculate data likelihood
            for (int k=0; k < maxK; ++k) {
                loglh = 0;
                if (debug) Rcout << "Calculating likelihood\n";
                if (debug) Rcout << "k = " << k << "\n";
                for (int d=0; d < P; ++d) {
                    loglh += df(i, d) * log(thisTheta(k, d)) + (1 - df(i, d)) * log(1 - thisTheta(k, d));
                    if (debug) {
                        Rcout << "d = " << d << "\n";
                        Rcout << "Theta val: " << thisTheta(k, d) << "\tx: " << df(i, d) << "\n";
                        Rcout << "Updated loglh " << loglh << "\n";
                    }
                }

                // Then calculate probabilities per cluster
                dummy = exp(log(pi_sampled(j-1, k)) + loglh);
                s[k] = dummy;
                cum_probs += dummy;
            }

            if (debug) {
                Rcout << "Raw probs:\n";
                for (int p = 0; p < maxK; ++p) {
                    Rcout << s[p] << " | ";
                }
                Rcout << "\n";
            }

            // Then normalise
            for (int p = 0; p < maxK; ++p) {
                s[p] /= cum_probs;
            }

            if (debug) {
                Rcout << "Normalised probs:\n";
                for (int p = 0; p < maxK; ++p) {
                    Rcout << s[p] << " | ";
                }
                Rcout << "\n";
            }

            // Now can draw labels
            rmultinom(1, s.begin(), maxK, this_z.begin());
            if (debug) Rcout << "this_z: " << this_z << "\n";
            z_sampled.slice(j).row(i) = as<arma::vec>(this_z).t();
        }  // End looping through individuals

        if (debug) Rcout << "\n\nNow going to sample pi";
        
        // Simulate V values
        
        // Now calculate number of patients in each cluster and sum of data points as before
        IntegerVector ck(maxK);
        IntegerMatrix Vkd(maxK, P);
        // Doing this in reverse as need number in previous clusters for 
        // stick breaking
        int num_previous_clusters = 0;
        arma::vec v(maxK);
        for (int k = maxK-1; k >= 0; --k) {
            // Get number of people in each cluster
            for (int i=0; i < N; ++i) {
                Znk = z_sampled(i, k, j);
                ck[k] += Znk;
                for (int d = 0; d < P; ++d) {
                    if (debug) {
                        Rcout << "k: " << k << "\t";
                        Rcout << "d: " << d << "\t";
                        Rcout << "i: " << i << "\t";
                        Rcout << "Znk: " << Znk << "\t";
                        Rcout << "ck[k]: " << ck[k] << "\t";
                        Rcout << "df[i, d]: " << df(i, d) << "\n";
                    }
                    Vkd(k, d) += Znk * df(i, d);
                }
            }
            double beta1, beta2;
            beta1 = 1 + ck[k];
            beta2 = alpha + num_previous_clusters;
            v(k) = R::rbeta(beta1, beta2);
            if (debug) Rcout << "beta1: " << beta1 << "\tbeta2: " << beta2 << "\tv: " << v(k) << "\n";
            
            num_previous_clusters += ck[k];
        }
        v(maxK-1) = 1;
        
        // Stick breaking to get pis
        this_pi(0) = v(0);
        double cumprod = 1 - v(0);
        for (int k = 1; k < maxK; k++) {
            this_pi(k) = cumprod * v(k);
            cumprod *= (1 - v(k));
        }
        if (debug) Rcout << "this_pi: " << this_pi << "\n";
        pi_sampled(j, _) = this_pi;

        // Generate theta(t)kd from Beta(γkd+vkd,δkd+uk−vkd)(for allk,d)
        for (int k = 0; k < maxK; ++k) {
            for (int d = 0; d < P; ++d) {
                if (debug) {
                    Rcout << "k: " << k << "\td: " << d << "\t";
                    Rcout << "Vkd: " << Vkd(k, d) << "\tck: " << ck[k] << "\t";
                }
                theta_row(k, d) = R::rbeta(beta + Vkd(k, d), gamma + ck[k] - Vkd(k, d));
            }
        }
        if (debug) Rcout << "Theta: " << theta_row << "\n";
        theta_sampled.slice(j) = as<arma::mat>(theta_row);
    }

    // Determine cluster labels for sampled values, as currently are in binary format
    for (int j = 0; j < nsamples; ++j) {
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < maxK; ++k) {
                if (z_sampled(i, k, j) == 1) {
                    z_out(j, i) = k + 1;
                    continue;
                }
            }
        }
    }
    
    List ret;
    ret["pi"] = pi_sampled(Range(burnin, nsamples-1), _);
    ret["z"] = z_out.tail_rows(nsamples-burnin);
    arma::cube thetas_post = theta_sampled.tail_slices(nsamples - burnin);
    ret["theta"] = thetas_post;
    return(ret);
}

