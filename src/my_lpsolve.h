// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

void lp_transbig_edit(int, int, double *, double *);
arma::Mat<int> my_lpsolve(arma::mat x);