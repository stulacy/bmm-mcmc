// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "my_lpsolve.h"

arma::Mat<int> my_stephens_batch(arma::cube p, bool debug);