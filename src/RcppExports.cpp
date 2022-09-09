// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_cross
arma::mat calc_cross(arma::mat cross, arma::mat within);
RcppExport SEXP _GPCERF_calc_cross(SEXP crossSEXP, SEXP withinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cross(crossSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type within(withinSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_cross(cross, within));
    return rcpp_result_gen;
END_RCPP
}
// arma_mm
arma::vec arma_mm(const arma::mat& m, const arma::vec& v);
RcppExport SEXP _GPCERF_arma_mm(SEXP mSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_mm(m, v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GPCERF_calc_cross", (DL_FUNC) &_GPCERF_calc_cross, 2},
    {"_GPCERF_arma_mm", (DL_FUNC) &_GPCERF_arma_mm, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPCERF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
