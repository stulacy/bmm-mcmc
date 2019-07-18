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

    arma::rowvec foo = arma::linspace<arma::rowvec>(0, 1, 2);
    arma::mat bar = arma::randn<arma::mat>(5, 2);
    arma::mat cat = arma::randn<arma::mat>(5, 2);
    Rcout << "foo: " << foo << "\tsize: " << arma::size(foo) << "\n\n";
    Rcout << "bar: " << bar << "\tsize: " << arma::size(bar) << "\n\n";
    Rcout << "cat: " << cat << "\tsize: " << arma::size(cat) << "\n\n";

    arma::mat inner_1(5, 2);
    for (int d=0; d<5; ++d) {
      for (int k=0; k < 2; ++k) {
        inner_1(d, k) = foo(0, k) - bar(d, k);
      }
    }


    Rcout << "inner1: " << inner_1 << "\tsize: " << arma::size(inner_1) << "\n";

    arma::mat inner_2(5, 2);
    for (int k=0; k < 2; ++k) {
      inner_2.col(k) =  - bar.col(k) + foo(0, k);
    }
    Rcout << "inner2: " << inner_2 << "\tsize: " << arma::size(inner_2) << "\n";

    Rcout << "bar + cat: " << bar+cat << "\n";

    //Rcout << "bar[foo]: " << bar.elem(foo) << "\tsize: " << arma::size(bar.elem(foo)) << "\n";

    //foo.shed_row(3);
    //Rcout << "foo: " << foo << "\tsize: " << arma::size(foo) << "\n";

  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
