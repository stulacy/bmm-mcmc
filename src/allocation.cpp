// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"

using namespace Rcpp;

// Allocation sampler from Nobile and Fearnside 2007
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.144.8651&rep=rep1&type=pdf
// [[Rcpp::export]]
List allocation_cpp(IntegerMatrix df,
                    IntegerVector initialK,
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

    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;

    // Random sample for first row
    arma::Mat<int> z_out(nsamples, N);
    arma::Mat<int> z_out_relabelled(nsamples, N);
    z_out.row(0) = as<arma::Row<int>>(initialK);
    arma::cube thetas(maxK, P, nsamples);
    arma::cube thetas_relab(maxK, P, nsamples);
    arma::vec alpha_sampled(nsamples);
    if (alpha == 0) {
        alpha_sampled(0) = 1;
    } else {
        alpha_sampled.fill(alpha);
    }

    int curr_cluster;

    // Vector for each cluster to track members
    // Helps to quickly identify Ck
    std::vector< std::vector< int > > clusters(maxK);
    for (int i=0; i < N; ++i) {
        clusters[initialK(i)-1].push_back(i);
    }
    std::vector< std::vector< int > > temp_clusters;

    std::vector <int> Ck;

    NumericVector probs(maxK);
    IntegerVector this_z(maxK);
    double probs_sum;
    int Nk, sum_xd, xnd, dsum;
    double LHS, logLH, left, right, denom, full, dummy;
    
    // Data structures for relabelling
    arma::cube probs_out(N, maxK, burnrelabel, arma::fill::zeros);
    arma::mat probs_sample(N, maxK, arma::fill::zeros);
    arma::mat Q;
    std::pair<arma::Row<int>, arma::mat> stephens_out;
    arma::Row<int> permutations_sample(maxK);
    arma::Mat<int> permutations(nsamples - burnin, maxK);
    
    int n_moves = 4;
    int move_choice;
    int K = 2;

    // At each sample, for each person:
    for (int j=1; j < nsamples; ++j) {
        
        // TODO Implement prior on K
        
        // Gibbs update
        if (debug) Rcout << "TODO: Implement Gibbs update\n";
        
        // Not implementing M3 as assume symmetric TODO add symmetry check
        //move_choice = arma::randperm(n_moves, 1)(0);
        move_choice = 1;
        switch (move_choice) {
            case 0: {
                Rcout << "M1\n";
                if (K < 2) {
                    break;
                }
                // Select 2 clusters
                arma::uvec to_perm = arma::randperm(K, 2);
                Rcout << "Selected clusters " << to_perm << "\n";
                
                // Calculate p1
                double p1 = R::rbeta(alpha_sampled(j-1), alpha_sampled(j-1));
                
                // Work out moves
                arma::Row<int> tempz = z_out.row(j-1);
                temp_clusters = clusters;
                
                // Find all inds in first clust
                for (int i : clusters[to_perm(0)]) {
                    if (arma::randu() < (1 - p1)) {
                        tempz(i) = to_perm(1) + 1;
                        temp_clusters[to_perm(0)].erase(std::remove(temp_clusters[to_perm(0)].begin(),
                                                                    temp_clusters[to_perm(0)].end(),
                                                                    i),
                                                        temp_clusters[to_perm(0)].end());
                        temp_clusters[to_perm(1)].push_back(i);
                    } 
                }
                for (int i : clusters[to_perm(1)]) {
                    if (arma::randu() < p1) {
                        tempz(i) = to_perm(0) + 1;
                        temp_clusters[to_perm(1)].erase(std::remove(temp_clusters[to_perm(1)].begin(),
                                                                    temp_clusters[to_perm(1)].end(),
                                                                    i),
                                                        temp_clusters[to_perm(1)].end());
                        temp_clusters[to_perm(0)].push_back(i);
                    }
                }
                
                // Calculate acceptance prob from Eq 6 in Nobile & Fearnside
                // which uses predictive distribution, see formulation here
                // using gamma function
                // https://en.wikipedia.org/wiki/Beta-binomial_distribution
                arma::uvec nNew(2);
                arma::uvec nOld(2);
                arma::mat sOld(2, P, arma::fill::zeros);
                arma::mat sNew(2, P, arma::fill::zeros);
                for (int i = 0; i < 2; ++i) {
                    nOld(i) = clusters[to_perm(i)].size();
                    nNew(i) = temp_clusters[to_perm(i)].size();
                    for (int d = 0; d < P; ++d) {
                        for (int c : clusters[to_perm(i)]) {
                            sOld(i, d) += df_arma(c, d);
                        }
                        for (int c : temp_clusters[to_perm(i)]) {
                            sNew(i, d) += df_arma(c, d);
                        }
                    }
                }
                double logAR = 0;
                for (int i = 0; i < 2; ++i) {
                    logAR += arma::sum(lgamma(beta + sNew.row(i)) + lgamma(gamma + nNew(i) - sNew.row(i)) - lgamma(beta + sOld.row(i)) - lgamma(gamma + nOld(i) - sOld.row(i))) + P * (lgamma(beta + gamma + nOld(i)) - lgamma(beta + gamma + nNew(i)));
                }
                Rcout << "logAR: " << logAR << "\n";
                if (log(arma::randu()) < logAR) {
                    Rcout << "Accepted\n";
                    clusters = temp_clusters;
                    z_out = tempz;
                }
                break;
            }
            case 1: {
                Rcout << "M2\n";
                if (K < 2) {
                    break;
                }
                
                // Select 2 clusters
                arma::uvec to_perm = arma::randperm(K, 2);
                
                // Select number of inds m from uniform {1, n_j1}
                int ninds = arma::randperm(clusters[to_perm(0)].size(), 1)[0] + 1;
                arma::uvec inds_perm = arma::randperm(clusters[to_perm(0)].size(), ninds);
                
                Rcout << "Chosen m individuals: " << ninds << "\n";
                
                // Move these individuals to component j2
                arma::Row<int> tempz = z_out.row(j-1);
                temp_clusters = clusters;
                
                // Find all inds in first clust
                for (int i : inds_perm) {
                    tempz(i) = to_perm(1) + 1;
                    temp_clusters[to_perm(0)].erase(std::remove(temp_clusters[to_perm(0)].begin(),
                                                                temp_clusters[to_perm(0)].end(),
                                                                i),
                                                    temp_clusters[to_perm(0)].end());
                    temp_clusters[to_perm(1)].push_back(i);
                }
                
                // Calculate s values for acceptance ratio
                arma::uvec nNew(2);
                arma::uvec nOld(2);
                arma::mat sOld(2, P, arma::fill::zeros);
                arma::mat sNew(2, P, arma::fill::zeros);
                for (int i = 0; i < 2; ++i) {
                    nOld(i) = clusters[to_perm(i)].size();
                    nNew(i) = temp_clusters[to_perm(i)].size();
                    for (int d = 0; d < P; ++d) {
                        for (int c : clusters[to_perm(i)]) {
                            sOld(i, d) += df_arma(c, d);
                        }
                        for (int c : temp_clusters[to_perm(i)]) {
                            sNew(i, d) += df_arma(c, d);
                        }
                    }
                }
                double logAR = 0;
                break;
            }
            case 2: {
                Rcout << "M3\n";
                break;
            }
            case 3: {
                Rcout << "M4\n";
                break;
            }
            default: {
                Rcerr << "Uh-oh, shouldn't end up here\n";
            }
        }
        
    }
    
    // Due to how subcubes work (i.e. when have sliced them), can't just
    // return directly, but instead need to force as cube
    List ret;
    //arma::cube thetas_post = thetas.tail_slices(nsamples - burnin);
    //ret["alpha"] = alpha_sampled.tail_rows(nsamples-burnin);
    //ret["permutations"] = permutations;
    //if (relabel) {
    //    arma::cube thetas_relabelled = thetas_relab.tail_slices(nsamples - burnin);
    //    ret["z"] = z_out_relabelled.tail_rows(nsamples-burnin);
    //    ret["theta"] = thetas_relabelled;
    //    ret["z_original"] = z_out.tail_rows(nsamples-burnin);
    //    ret["theta_original"] = thetas_post;
    //} else {
    //    ret["z"] = z_out.tail_rows(nsamples-burnin);
    //    ret["theta"] = thetas_post;
    //}
    return ret;
}

