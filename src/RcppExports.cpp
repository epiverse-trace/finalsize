// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>& Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// final_size_cpp
Eigen::VectorXd final_size_cpp(const double& r0,
                               const Eigen::MatrixXd& contact_matrix,
                               const Eigen::VectorXd& demography_vector,
                               const Eigen::VectorXd& prop_initial_infected,
                               Eigen::VectorXd prop_suscep,
                               const int iterations);
RcppExport SEXP _finalsize_final_size_cpp(SEXP r0SEXP, SEXP contact_matrixSEXP,
                                          SEXP demography_vectorSEXP,
                                          SEXP prop_initial_infectedSEXP,
                                          SEXP prop_suscepSEXP,
                                          SEXP iterationsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const double&>::type r0(r0SEXP);
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type contact_matrix(
      contact_matrixSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type demography_vector(
      demography_vectorSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type
      prop_initial_infected(prop_initial_infectedSEXP);
  Rcpp::traits::input_parameter<Eigen::VectorXd>::type prop_suscep(
      prop_suscepSEXP);
  Rcpp::traits::input_parameter<const int>::type iterations(iterationsSEXP);
  rcpp_result_gen = Rcpp::wrap(
      final_size_cpp(r0, contact_matrix, demography_vector,
                     prop_initial_infected, prop_suscep, iterations));
  return rcpp_result_gen;
  END_RCPP
}
// final_size_grps_cpp
Eigen::ArrayXd final_size_grps_cpp(const Eigen::MatrixXd& contact_matrix,
                                   const Eigen::VectorXd& demography_vector,
                                   const Eigen::MatrixXd& p_susceptibility,
                                   const Eigen::MatrixXd& susceptibility);
RcppExport SEXP _finalsize_final_size_grps_cpp(SEXP contact_matrixSEXP,
                                               SEXP demography_vectorSEXP,
                                               SEXP p_susceptibilitySEXP,
                                               SEXP susceptibilitySEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type contact_matrix(
      contact_matrixSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type demography_vector(
      demography_vectorSEXP);
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type p_susceptibility(
      p_susceptibilitySEXP);
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type susceptibility(
      susceptibilitySEXP);
  rcpp_result_gen = Rcpp::wrap(final_size_grps_cpp(
      contact_matrix, demography_vector, p_susceptibility, susceptibility));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_finalsize_final_size_cpp", (DL_FUNC)&_finalsize_final_size_cpp, 6},
    {"_finalsize_final_size_grps_cpp", (DL_FUNC)&_finalsize_final_size_grps_cpp,
     4},
    {"run_testthat_tests", (DL_FUNC)&run_testthat_tests, 1},
    {NULL, NULL, 0}};

RcppExport void R_init_finalsize(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
