// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "my_lpsolve.h"

arma::mat my_stephens_batch(arma::cube, bool);
std::pair<arma::Row<int>, arma::mat> my_stephens_online(arma::mat, arma::mat, int, bool);