// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "my_lpsolve.h"

arma::mat my_stephens_batch(arma::cube p, bool debug);