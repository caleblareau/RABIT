// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hellow
List rcpp_hellow(int q);
RcppExport SEXP RABIT_rcpp_hellow(SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    __result = Rcpp::wrap(rcpp_hellow(q));
    return __result;
END_RCPP
}
