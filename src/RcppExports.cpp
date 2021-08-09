// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// convert_methy_to_dss_cpp
std::vector<std::string> convert_methy_to_dss_cpp(std::string input, std::string output_dir);
RcppExport SEXP _NanoMethViz_convert_methy_to_dss_cpp(SEXP inputSEXP, SEXP output_dirSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_dir(output_dirSEXP);
    rcpp_result_gen = Rcpp::wrap(convert_methy_to_dss_cpp(input, output_dir));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NanoMethViz_convert_methy_to_dss_cpp", (DL_FUNC) &_NanoMethViz_convert_methy_to_dss_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_NanoMethViz(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
