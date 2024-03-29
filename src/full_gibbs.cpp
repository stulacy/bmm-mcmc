// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"
#include "stephens.h"

using namespace Rcpp;

// Taken from https://www.mjdenny.com/blog.html
// draws from dirichlet using relationship with gamma(alpha, 1)
// [[Rcpp::export]]
arma::vec rdirichlet_cpp(arma::vec alpha_m) {
    int distribution_size = alpha_m.n_elem;
    // each row will be a draw from a Dirichlet
    arma::vec distribution = arma::zeros(distribution_size);

    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
        double cur = R::rgamma(alpha_m[j], 1.0);
        distribution(j) = cur;
        sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
        distribution(j) = distribution(j) / sum_term;
    }
    return(distribution);
}

// Using Algorithm 2 from
// https://www4.stat.ncsu.edu/~wilson/prelim/Review1.pdf
// [[Rcpp::export]]
List gibbs_cpp(IntegerMatrix df,
                        NumericVector initialPi,
                        NumericMatrix initialTheta,
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

    int N = df.nrow();
    int P = df.ncol();


    // Random sample for first row
    NumericMatrix pi_sampled(nsamples, K);
    arma::cube z_sampled(N, K, nsamples);
    arma::cube theta_sampled(K, P, nsamples);
    arma::cube theta_relab(K, P, nsamples);
    arma::Mat<int> z_out(nsamples, N);
    arma::Mat<int> z_out_relabelled(nsamples, N);
    
    // Structures for relabelling
    arma::cube probs_out(N, K, burnrelabel, arma::fill::zeros);
    arma::mat probs_sample(N, K, arma::fill::zeros);
    arma::mat Q;
    std::pair<arma::Row<int>, arma::mat> stephens_out;
    arma::Row<int> permutations_sample(K);
    arma::Mat<int> permutations(nsamples - burnin, K);

    int Znk;
    double loglh, cum_probs, dummy, theta_update;
    pi_sampled(0, _) = initialPi;
    theta_sampled.slice(0) = as<arma::mat>(initialTheta);
    NumericMatrix thisTheta(K,P);
    IntegerVector this_z(K);
    NumericVector s(K);
    arma::vec dirich_params = arma::zeros(K);
    NumericVector this_pi(K);
    arma::vec alpha_sampled(nsamples);
    if (alpha == 0) {
        alpha_sampled(0) = 1;
    } else {
        alpha_sampled.fill(alpha);
    }

    for (int j=1; j < nsamples; ++j) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";

        thisTheta = wrap(theta_sampled.slice(j-1));
        for (int i=0; i < N; ++i) {
            if (debug) Rcout << "Iterating through N to sample z\n";
            cum_probs = 0;

            // Calculate data likelihood
            for (int k=0; k < K; ++k) {
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
                for (int p = 0; p < K; ++p) {
                    Rcout << s[p] << " | ";
                }
                Rcout << "\n";
            }

            // Then normalise
            for (int p = 0; p < K; ++p) {
                s[p] /= cum_probs;
            }

            if (debug) {
                Rcout << "Normalised probs:\n";
                for (int p = 0; p < K; ++p) {
                    Rcout << s[p] << " | ";
                }
                Rcout << "\n";
            }

            // Now can draw labels
            rmultinom(1, s.begin(), K, this_z.begin());
            if (debug) Rcout << "this_z: " << this_z << "\n";
            z_sampled.slice(j).row(i) = as<arma::vec>(this_z).t();
            // Determine cluster labels for sampled values, as currently are in binary format
            for (int k = 0; k < K; ++k) {
                if (this_z(k) == 1) {
                    z_out(j, i) = k + 1;
                    continue;
                }
            }
            
            // Save initial probabilities so can do batch Stephens
            // to generate initial Q
            if (relabel) {
                if (j < burnin && j >= (burnin - burnrelabel)) {
                    for (int k=0; k < K; ++k) {
                        probs_out(i, k, j - burnin + burnrelabel) = s(k);
                    }
                } else if (j >= burnin) {
                    for (int k=0; k < K; ++k) {
                        probs_sample(i, k) = s(k);
                    }
                }
            }
        }
        
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
        

        if (debug) Rcout << "\n\nNow going to sample pi and thetas\n";
        // Now sample pi and thetas
        // Now calculate number of patients in each cluster and sum of data points as before
        IntegerVector ck(K);
        IntegerMatrix Vkd(K, P);
        for (int k = 0; k < K; ++k) {
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
        }

        for (int k = 0; k < K; k++) {
            dirich_params(k) = (alpha_sampled(j-1) / K) + ck[k];
        }
        if (debug) Rcout << "dirich params: " << dirich_params << "\n";

        // Generate pi(t) from Dirichlet(α1+u1,...,αK+uK)
        this_pi = wrap(rdirichlet_cpp(dirich_params));
        if (debug) Rcout << "this_pi: " << this_pi << "\n";
        pi_sampled(j, _) = this_pi;

        // Generate theta(t)kd from Beta(γkd+vkd,δkd+uk−vkd)(for allk,d)
        for (int k = 0; k < K; ++k) {
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
        
        // Update alpha
        if (alpha == 0) {
            alpha_sampled(j) = update_alpha(alpha_sampled(j-1), a, b, N, K);
        }
    }

    List ret;
    arma::cube thetas_post = theta_sampled.tail_slices(nsamples - burnin);
    ret["pi"] = pi_sampled(Range(burnin, nsamples-1), _);
    ret["alpha"] = alpha_sampled.tail_rows(nsamples-burnin);
    ret["permutations"] = permutations;
    if (relabel) {
        arma::cube thetas_relabelled = theta_relab.tail_slices(nsamples - burnin);
        ret["z"] = z_out_relabelled.tail_rows(nsamples-burnin);
        ret["theta"] = thetas_relabelled;
        ret["z_original"] = z_out.tail_rows(nsamples-burnin);
        ret["theta_original"] = thetas_post;
    } else {
        ret["z"] = z_out.tail_rows(nsamples-burnin);
        ret["theta"] = thetas_post;
    }
    return(ret);
}

