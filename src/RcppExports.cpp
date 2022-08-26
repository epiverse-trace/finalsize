// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// solve_final_size_internal
Eigen::MatrixXd solve_final_size_internal(const Eigen::MatrixXd& contact_matrix, const Eigen::VectorXd& demography, const Eigen::MatrixXd& p_susceptibility, const Eigen::MatrixXd& susceptibility, const bool adapt_step, const double tolerance);
RcppExport SEXP _finalsize_solve_final_size_internal(SEXP contact_matrixSEXP, SEXP demographySEXP, SEXP p_susceptibilitySEXP, SEXP susceptibilitySEXP, SEXP adapt_stepSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type contact_matrix(contact_matrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type demography(demographySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type p_susceptibility(p_susceptibilitySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type susceptibility(susceptibilitySEXP);
    Rcpp::traits::input_parameter< const bool >::type adapt_step(adapt_stepSEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_final_size_internal(contact_matrix, demography, p_susceptibility, susceptibility, adapt_step, tolerance));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_finalsize_solve_final_size_internal", (DL_FUNC) &_finalsize_solve_final_size_internal, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_finalsize(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
