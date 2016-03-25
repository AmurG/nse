// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// f_bootstrap
Rcpp::NumericVector f_bootstrap(Rcpp::NumericVector x, double p, int type);
RcppExport SEXP nse_f_bootstrap(SEXP xSEXP, SEXP pSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    __result = Rcpp::wrap(f_bootstrap(x, p, type));
    return __result;
END_RCPP
}