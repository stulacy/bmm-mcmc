// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// collapsed_gibbs_cpp
List collapsed_gibbs_cpp(IntegerMatrix df, IntegerVector initialK, int nsamples, int K, double alpha, double beta, double gamma, int burnin, bool debug);
RcppExport SEXP _bmmmcmc_collapsed_gibbs_cpp(SEXP dfSEXP, SEXP initialKSEXP, SEXP nsamplesSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP burninSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type df(dfSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type initialK(initialKSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_gibbs_cpp(df, initialK, nsamples, K, alpha, beta, gamma, burnin, debug));
    return rcpp_result_gen;
END_RCPP
}
// collapsed_gibbs_dp_cpp
List collapsed_gibbs_dp_cpp(IntegerMatrix df, int nsamples, double alpha, double beta, double gamma, double a, double b, int burnin, bool debug);
RcppExport SEXP _bmmmcmc_collapsed_gibbs_dp_cpp(SEXP dfSEXP, SEXP nsamplesSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP burninSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(collapsed_gibbs_dp_cpp(df, nsamples, alpha, beta, gamma, a, b, burnin, debug));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_cpp
arma::vec rdirichlet_cpp(arma::vec alpha_m);
RcppExport SEXP _bmmmcmc_rdirichlet_cpp(SEXP alpha_mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha_m(alpha_mSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(alpha_m));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_cpp
List gibbs_cpp(IntegerMatrix df, NumericVector initialPi, NumericMatrix initialTheta, int nsamples, int K, double alpha, double beta, double gamma, int burnin, bool debug);
RcppExport SEXP _bmmmcmc_gibbs_cpp(SEXP dfSEXP, SEXP initialPiSEXP, SEXP initialThetaSEXP, SEXP nsamplesSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP burninSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initialPi(initialPiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initialTheta(initialThetaSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_cpp(df, initialPi, initialTheta, nsamples, K, alpha, beta, gamma, burnin, debug));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _bmmmcmc_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bmmmcmc_collapsed_gibbs_cpp", (DL_FUNC) &_bmmmcmc_collapsed_gibbs_cpp, 9},
    {"_bmmmcmc_collapsed_gibbs_dp_cpp", (DL_FUNC) &_bmmmcmc_collapsed_gibbs_dp_cpp, 9},
    {"_bmmmcmc_rdirichlet_cpp", (DL_FUNC) &_bmmmcmc_rdirichlet_cpp, 1},
    {"_bmmmcmc_gibbs_cpp", (DL_FUNC) &_bmmmcmc_gibbs_cpp, 10},
    {"_bmmmcmc_timesTwo", (DL_FUNC) &_bmmmcmc_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_bmmmcmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
