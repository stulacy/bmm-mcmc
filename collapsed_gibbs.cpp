// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

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

// Collapsed Gibbs sampler
// For each i in 1:N
//  Sample z_i given (theta, z_(-i), x)
// For k in 1:K
//  Sample theta_k given (theta_(-k), z, x)
//
// Tu in Section 3.23 says needs to sample p(z | ...), p(pi | ...), and p(theta | ...)
// But if use collapsed Gibbs then can sample p(z) directly
// https://people.eecs.berkeley.edu/~stephentu/writeups/mixturemodels.pdf
//
// Eq 21 from Van Maarten also suggests can just use collapsed Gibbs sampler over z
// this equation looks to be the same as Tu's version.
// https://pdfs.semanticscholar.org/525c/ff658d34ae4d47e84b8ec4ede3ce6c561afc.pdf
//
// And this looks to be the same as 3.5 from Neal
// http://www.stat.columbia.edu/npbayes/papers/neal_sampling.pdf
// [[Rcpp::export]]
List collapsed_gibbs_cpp(IntegerMatrix df,
                         IntegerMatrix initialK,
                         int nsamples,
                         int K,
                         double alpha,
                         double beta,
                         double gamma,
                         bool debug) {

    arma::Mat<int> df_arma = as<arma::Mat<int>>(df);

    int N = df_arma.n_rows;
    int P = df_arma.n_cols;

    // Random sample for first row
    arma::Cube<int> allocations(N, K, nsamples);
    allocations.slice(0) = as<arma::Mat<int>>(initialK);
    arma::cube thetas(K, P, nsamples);

    arma::uvec ck_inds(N);
    arma::Mat<int> Ck_(N-1, K);
    arma::Mat<int> this_df(N-1, P);
    arma::vec N_nk(K);
    arma::mat sum_k_d(P, K);
    arma::mat inner1(P, K);
    arma::mat inner2(P, K);
    arma::vec left_term(K);
    arma::mat numerator(P, K);
    arma::vec denominator(K);
    arma::vec frac(K);
    arma::vec probs(K);
    arma::mat this_theta(K, P);
    arma::vec Nk(K);

    //arma::colvec N_nk;
    for (int j=1; j < nsamples; ++j) {
        Rcpp::Rcout << "Sample: " << j+1 << "\n";
        // Set params
        for (int i=0; i < N; ++i) {
            if (debug) Rcpp::Rcout << "\nObservation: " << i << "\n~~~~~~~~~~~~~~\n\n";
            ck_inds = arma::linspace<arma::uvec>(0, N-1, N);
            ck_inds.shed_row(i);
            Ck_ = allocations.slice(j-1).rows(ck_inds);
            this_df = df_arma.rows(ck_inds);
            N_nk = arma::conv_to<arma::vec>::from(arma::sum(Ck_, 0));
            if (debug) {
                Rcout << "Ck_ " << Ck_ << "\n";
                Rcout << "Ck_ size: " << arma::size(Ck_) << "\n";
                Rcout << "N_nk " << N_nk << "\n";
                Rcout << "N_nk size " << arma::size(N_nk) << "\n";
                Rcout << "df size: " << arma::size(this_df) << "\n";
            }

            sum_k_d = arma::conv_to<arma::mat>::from(this_df.t() * Ck_);
            if (debug) {
                Rcout << "sum_k_d: " << sum_k_d << "\n";
                Rcout << "size: sum_k_d " << arma::size(sum_k_d) << "\n";
            }

            // ***** LEFT HAND SIDE of density *****
            left_term = arma::log((N_nk + (alpha / K)) / (N - 1 + alpha));
            if (debug) Rcout << "left_term: " << left_term << "\tsize: " << arma::size(left_term) << "\n";

            // ***** RIGHT HAND SIDE of density *****
            // Calculate inner parts of numerator, both returning D x K matrices
            inner1 = arma::log(beta + sum_k_d);

            // Now transform sum_k_d so can use it for right hand inner sum
            sum_k_d.each_row() -= N_nk.t();
            inner2 = arma::log(gamma - sum_k_d);
            if (debug) Rcout << "sum_k_d after transform: " << sum_k_d << "\tsize: " << arma::size(sum_k_d) << "\n";
            if (debug) Rcout << "inner1: " << inner1 << "\tsize: " << arma::size(inner1) << "\n";
            if (debug) Rcout << "inner2: " << inner2 << "\tsize: " << arma::size(inner2) << "\n";

            // Add inner sides up and divide by denominator to get DxK for this individual
            inner1.each_col() %= arma::conv_to<arma::vec>::from(df_arma.row(i).t());
            inner2.each_col() %= arma::conv_to<arma::vec>::from(1-df_arma.row(i).t());
            if (debug) Rcout << "inner1 multiplied by data: " << inner1 << "\tsize: " << arma::size(inner1) << "\n";
            if (debug) Rcout << "inner2 multiplied by data: " << inner2 << "\tsize: " << arma::size(inner2) << "\n";

            // TODO DO THIS RHS IN ONE STEP
            numerator = inner1 + inner2;
            if (debug) Rcout << "numerator: " << numerator << "\tsize: " << arma::size(numerator) << "\n";

            // Denominator is also K length vector due to N_nk
            denominator = arma::log(beta + gamma + N_nk);
            if (debug) Rcout << "denominator: " << denominator << "\tsize: " << arma::size(denominator) << "\n";

            // Element wise division to obtain K length vector for right-hand side of density
            numerator.each_row() -= denominator.t();
            if (debug) Rcout << "numerator: " << numerator << "\tsize: " << arma::size(numerator) << "\n";

            // Sum over D and add to left term to get raw probabilities for this observation
            probs = arma::exp(left_term + arma::sum(numerator, 0).t());
            if (debug) Rcout << "Raw probs: " << probs << "\tsize: " << arma::size(probs) << "\n";
            probs /= arma::sum(probs);
            if (debug) Rcout << "normalised probs: " << probs << "\tsize: " << arma::size(probs) << "\n";

            // Now pull out the Z assignments
            IntegerVector this_z(K);
            R::rmultinom(1, probs.memptr(), K, this_z.begin());
            if (debug) Rcout << "this_z: " << this_z << "\n";
            allocations.slice(j).row(i) = as<arma::Row<int>>(this_z);
        }
        // Calculate thetas here to be returned
        this_theta = (arma::conv_to<arma::mat>::from(df_arma.t() * allocations.slice(j))).t();
        if (debug) Rcout << "this_theta: " << this_theta << "\tsize: " << arma::size(this_theta) << "\n";
        Nk = arma::conv_to<arma::vec>::from(arma::sum(allocations.slice(j), 0));
        if (debug) Rcout << "Nk: " << Nk << "\tsize: " << arma::size(Nk) << "\n";

        this_theta.each_col() /= Nk;
        this_theta.replace(arma::datum::nan, 0);
        if (debug) Rcout << "this_theta: " << this_theta << "\tsize: " << arma::size(this_theta) << "\n";

        thetas.slice(j) = this_theta;
    }

    List ret;
    ret["z"] = allocations;
    ret["theta"] = thetas;
    return ret;
}

