// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"
#include "stephens.h"

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
                             double a,
                             double b,
                             int burnin,
                             bool relabel,
                             int burnrelabel,
                             bool debug) {

    int N = df.nrow();
    int P = df.ncol();


    // Random sample for first row
    NumericMatrix pi_sampled(nsamples, maxK);
    arma::cube z_sampled(N, maxK, nsamples);
    arma::cube theta_sampled(maxK, P, nsamples);
    arma::cube theta_relab(maxK, P, nsamples);
    arma::Mat<int> z_out(nsamples, N);
    arma::Mat<int> z_out_relabelled(nsamples, N);

    int Znk;
    // If pi is greater than this threshold then we deem it a 'viable'
    // cluster. Number of viable clusters only used for updating alpha
    double viable_threshold = 0.01;
    int K_viable = maxK;
    double loglh, cum_probs, dummy, theta_update;
    pi_sampled(0, _) = initialPi;
    arma::vec alpha_sampled(nsamples);
    alpha_sampled(0) = alpha;
    theta_sampled.slice(0) = as<arma::mat>(initialTheta);
    NumericMatrix thisTheta(maxK,P);
    IntegerVector this_z(maxK);
    NumericVector s(maxK);
    arma::vec dirich_params = arma::zeros(maxK);
    NumericVector this_pi(maxK);
    NumericMatrix theta_row(maxK, P);
    arma::cube probs_out(N, maxK, burnrelabel, arma::fill::zeros);
    arma::mat probs_sample(N, maxK, arma::fill::zeros);
    arma::mat Q;
    std::pair<arma::Row<int>, arma::mat> stephens_out;
    arma::Row<int> permutations_sample(maxK);
    arma::Mat<int> permutations(nsamples - burnin, maxK);

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
            // Determine cluster labels for sampled values, as currently are in binary format
            for (int k = 0; k < maxK; ++k) {
                if (this_z(k) == 1) {
                    z_out(j, i) = k + 1;
                    continue;
                }
            }
            
            // Save initial probabilities so can do batch Stephens
            // to generate initial Q
            if (relabel) {
                if (j < burnin && j >= (burnin - burnrelabel)) {
                    for (int k=0; k < maxK; ++k) {
                        probs_out(i, k, j - burnin + burnrelabel) = s(k);
                    }
                } else if (j >= burnin) {
                    for (int k=0; k < maxK; ++k) {
                        probs_sample(i, k) = s(k);
                    }
                }
            }
        }  // End looping through individuals
        
        // To relabel clusters use Stephen's 2000b online algorithm.
        // Firstly need to initialise Q with values taken from a batch formulation
        // over $burnrelabel samples
        if (relabel) {
            if (j == (burnin - 1)) {
                Rcout << "Running Stephen's batch relabelling algorithm.\n";
                Q = my_stephens_batch(probs_out, false);
            } else if (j >= burnin) {
                stephens_out = my_stephens_online(Q, probs_sample, j, false);
                Q = stephens_out.second;
                permutations_sample = stephens_out.first;
                // Relabel Z 
                permutations.row(j-burnin) = permutations_sample;
                for (int i = 0; i < N; ++i) {
                    z_out_relabelled(j, i) = permutations_sample(z_out(j, i)-1) + 1;
                }
            }
        }
        

        if (debug) Rcout << "\n\nNow going to sample pi";
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
            beta2 = alpha_sampled(j-1) + num_previous_clusters;
            v(k) = R::rbeta(beta1, beta2);
            if (debug) Rcout << "beta1: " << beta1 << "\tbeta2: " << beta2 << "\tv: " << v(k) << "\n";
            
            num_previous_clusters += ck[k];
        }
        v(maxK-1) = 1;
        
        
        // TODO Can clean this up, move k=0 into loop
        K_viable = 0;
        // Stick breaking to get pis
        this_pi(0) = v(0);
        if (this_pi(0) > viable_threshold) {
            K_viable++;
        }
        double cumprod = 1 - v(0);
        for (int k = 1; k < maxK; k++) {
            this_pi(k) = cumprod * v(k);
            if (this_pi(k) > viable_threshold) {
                K_viable++;
            }
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
                theta_update = R::rbeta(beta + Vkd(k, d), gamma + ck[k] - Vkd(k, d));
                theta_sampled(k, d, j) = theta_update;
                if (relabel && j >= burnin) {
                    theta_relab(permutations_sample(k), d, j) = theta_update;
                }
            }
        }
        if (debug) Rcout << "Theta: " << theta_row << "\n";
        
        // Update alpha
        alpha_sampled(j) = update_alpha(alpha_sampled(j-1), a, b, N, K_viable);
    }  // End sampling loop
    
    List ret;
    arma::cube thetas_post = theta_sampled.tail_slices(nsamples - burnin);
    ret["pi"] = pi_sampled(Range(burnin, nsamples-1), _);
    ret["alpha"] = alpha_sampled.tail_rows(nsamples-burnin);
    if (relabel) {
        ret["z"] = z_out_relabelled.tail_rows(nsamples-burnin);
        arma::cube thetas_relabelled = theta_relab.tail_slices(nsamples - burnin);
        ret["theta"] = thetas_relabelled;
        ret["z_original"] = z_out.tail_rows(nsamples-burnin);
        ret["theta_original"] = thetas_post;
    } else {
        ret["z"] = z_out.tail_rows(nsamples-burnin);
        ret["theta"] = thetas_post;
    }
    
    return(ret);
}

