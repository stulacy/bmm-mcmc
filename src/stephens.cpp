#include "stephens.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat my_stephens_batch(arma::cube p, bool debug) {
  
    int N = p.n_rows;
    int K = p.n_cols;
    int M = p.n_slices;
    
    arma::mat cost_matrix(K, K);
    arma::Mat<int> perm(M, K);
    
    arma::Col<int> vec(M);
    for (int k=0; k<K; ++k) {
        vec.fill(k);
        perm.col(k) = vec;
    }
    
    double previous = -99;
    double current;
    double criterion = 99;
    double threshold = 10^(-6);
    int maxiter = 100;
    int t = 0;
    arma::mat temp(N, K);
    arma::Mat<int> solution(K, K);
  	arma::mat q(N, K);
    double min_prob = 0.000001;
    p.replace(0, min_prob);
    
    while((criterion > threshold) && (t < maxiter)) {
      	t++;
      	
        //compute q matrix
      	q.fill(0);
      	for (int k=0; k < K; ++k) {
      		for (int iter=0; iter < M; ++iter) {
      		    q.col(k) += p.slice(iter).col(perm(iter, k));
  		    }
      	}
      	q /= M;
      	
        arma::mat sub(N, K);
    	
		for (int iter=0; iter < M; ++iter) {
          	for (int k=0; k < K; ++k) {
          	    sub = arma::log(p.slice(iter));
          	    cost_matrix.row(k) = arma::sum(p.slice(iter) % (sub.each_col() - arma::log(q.col(k))), 0);
    		}
    		solution = my_lpsolve(cost_matrix);
    		for (int k=0; k < K; ++k) {
    		  perm(iter,k) = index_max(solution.col(k));
    		}
    		perm.row(iter) <- arma::sort_index(perm.row(iter));
    	}
    	current = arma::accu(cost_matrix % solution);
    	criterion = abs(previous - current);
    	if (debug) Rcout << "Iteration: " << t << "\tCost: " << criterion << "\tCurrent: " << current << "\tprevious: " << previous << "\n";
    	previous = current;
    }
    return q;
}

std::pair<arma::Row<int>, arma::mat> my_stephens_online(arma::mat q, arma::mat p, 
                                                        int sample_num, bool debug) {
  
    int N = p.n_rows;
    int K = p.n_cols;
    
    arma::mat cost_matrix(K, K);
    arma::Col<int> perm(K);
    arma::Mat<int> solution(K, K);
    arma::mat p_reordered(N, K);
    arma::mat q_new(N, K);
    
  	for (int k=0; k < K; ++k) {
  	    cost_matrix.row(k) = arma::sum(p % (p.each_col() - arma::log(q.col(k))), 0);
	}
	solution = my_lpsolve(cost_matrix);
	for (int k=0; k < K; ++k) {
	  perm(k) = index_max(solution.col(k));
	}
	perm <- arma::sort_index(perm);
  	if (debug) Rcout << "Perm: " << perm << "\n";
  	for (int k = 0; k < K; ++k) {
  	    p_reordered.col(k) = p.col(perm(k));
  	}
  	
  	// Finalise Q
  	q_new = (sample_num * (q + p_reordered)) / (sample_num + 1);
    return std::make_pair(perm.t(), q_new);
}