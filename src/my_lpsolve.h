// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdio.h>
#include "ls_source/lp_lib.h"

using namespace Rcpp;

void lp_transbig_edit(int, int, double *, double *);
arma::Mat<int> my_lpsolve(arma::mat x);