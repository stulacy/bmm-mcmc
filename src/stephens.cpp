#include "stephens.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::Mat<int> my_stephens_batch(arma::cube p, bool debug) {
  
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
    int maxiter = 8;
    int t = 0;
    arma::mat temp(N, K);
    arma::Mat<int> solution(K, K);
    
    while((criterion > threshold) && (t < maxiter)) {
        if (debug) Rcout << "Iteration " << t << "\n";
      	t++;
      	
        //compute q matrix
      	arma::mat q(N, K);
      	q.fill(0);
      	for (int k=0; k < K; ++k) {
      		for (int iter=0; iter < M; ++iter) {
      		    q.col(k) += p.slice(iter).col(perm(iter, k));
  		    }
      	}
      	if (debug) Rcout << "q before sum" << q << "\n";
      	q /= M;
      	if (debug) Rcout << "q after sum" << q << "\n";
      	
    arma::mat sub(N, K);
    	
		for (int iter=0; iter < M; ++iter) {
      	for (int k=0; k < K; ++k) {
      	    if (debug) Rcout << "p.slice" << p.slice(iter) << "\n";
      	    if (debug) Rcout << "log p.slice" << arma::log(p.slice(iter)) << "\n";
      	    if (debug) Rcout << "q.col" << arma::log(q.col(k)) << "\n";
      	    if (debug) Rcout << "SORT q.col" << arma::sort_index(arma::log(q.col(k))) << "\n";
      	    sub = arma::log(p.slice(iter));
      	    if (debug) Rcout << "raw: " << p.slice(iter) % (sub.each_col() - arma::log(q.col(k))) << "\n";
      	    if (debug) Rcout << "sum: " << arma::sum(p.slice(iter) % (sub.each_col() - arma::log(q.col(k))), 0) << "\n";
      	    cost_matrix.row(k) = arma::sum(p.slice(iter) % (sub.each_col() - arma::log(q.col(k))), 0);
    		}
      	if (debug) Rcout << "Cost matrix: " << cost_matrix << "\n";
    		solution = my_lpsolve(cost_matrix);
      	if (debug) Rcout << "Solution: " << solution << "\n";
    		for (int k=0; k < K; ++k) {
    		  perm(iter,k) = index_max(solution.col(k));
    		}
    		perm.row(iter) <- arma::sort_index(perm.row(iter));
    	}
    	arma::mat foo = cost_matrix % solution;
		  if (debug) Rcout << "foo: " << foo << "\n";
    	current = arma::accu(foo);
		  if (debug) Rcout << "Current: " << current << "\n";
    	criterion = abs(previous - current);
    	previous = current;
    }
    return perm;
}