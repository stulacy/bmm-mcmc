// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <stdlib.h>
#include "src/lp_lib.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
    arma::mat foo(2,2);
    foo(0, 0) = -1.267420;
    foo(0, 1) = -1.267420;
    foo(1, 0) = 2.074419;
    foo(1, 1) = 2.074419;
    Rcout << "Foo: " << foo  << "\n";
    return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
