// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
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
  
    arma::uvec foo = arma::linspace<arma::uvec>(0, 4, 5);
    arma::mat bar = arma::randn<arma::mat>(5, 2);
    Rcout << "foo: " << foo << "\tsize: " << arma::size(foo) << "\n";
    Rcout << "bar: " << bar << "\tsize: " << arma::size(bar) << "\n";
    Rcout << "bar[foo]: " << bar.elem(foo) << "\tsize: " << arma::size(bar.elem(foo)) << "\n";
    
    foo.shed_row(3);
    Rcout << "foo: " << foo << "\tsize: " << arma::size(foo) << "\n";
    
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
