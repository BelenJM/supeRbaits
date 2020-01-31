// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getBaits
Rcpp::DataFrame getBaits(std::string db_path, Rcpp::DataFrame df);
RcppExport SEXP _supeRbaits_getBaits(SEXP db_pathSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type db_path(db_pathSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(getBaits(db_path, df));
    return rcpp_result_gen;
END_RCPP
}
// getChromLengths
DataFrame getChromLengths(std::string path);
RcppExport SEXP _supeRbaits_getChromLengths(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(getChromLengths(path));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_supeRbaits_getBaits", (DL_FUNC) &_supeRbaits_getBaits, 2},
    {"_supeRbaits_getChromLengths", (DL_FUNC) &_supeRbaits_getChromLengths, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_supeRbaits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
