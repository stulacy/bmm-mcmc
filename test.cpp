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
    std::vector<double> probs(3);
    probs[0] = 0.1;
    probs[1] = 0.2;
    probs[2] = 0.7;

    NumericVector probs2(3);
    probs2[0] = 0.1;
    probs2[1] = 0.2;
    probs2[2] = 0.7;

    IntegerVector ans(3);

    rmultinom(1, probs.begin(), 3, ans.begin());
    Rcout << ans << "\n";
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
