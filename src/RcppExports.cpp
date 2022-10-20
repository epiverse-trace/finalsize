// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>& Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// epi_spread
Rcpp::List epi_spread(const Eigen::MatrixXd& contact_matrix,
                      const Eigen::VectorXd& demography_vector,
                      const Eigen::MatrixXd& p_susceptibility,
                      const Eigen::MatrixXd& susceptibility);
RcppExport SEXP _finalsize_epi_spread(SEXP contact_matrixSEXP,
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
  rcpp_result_gen = Rcpp::wrap(epi_spread(contact_matrix, demography_vector,
                                          p_susceptibility, susceptibility));
  return rcpp_result_gen;
  END_RCPP
}
// solve_final_size_iterative
Eigen::ArrayXd solve_final_size_iterative(
    const Eigen::MatrixXd& contact_matrix,
    const Eigen::VectorXd& demography_vector,
    const Eigen::VectorXd& susceptibility, const int iterations,
    const double tolerance, double step_rate, const bool adapt_step);
RcppExport SEXP _finalsize_solve_final_size_iterative(
    SEXP contact_matrixSEXP, SEXP demography_vectorSEXP,
    SEXP susceptibilitySEXP, SEXP iterationsSEXP, SEXP toleranceSEXP,
    SEXP step_rateSEXP, SEXP adapt_stepSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type contact_matrix(
      contact_matrixSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type demography_vector(
      demography_vectorSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type susceptibility(
      susceptibilitySEXP);
  Rcpp::traits::input_parameter<const int>::type iterations(iterationsSEXP);
  Rcpp::traits::input_parameter<const double>::type tolerance(toleranceSEXP);
  Rcpp::traits::input_parameter<double>::type step_rate(step_rateSEXP);
  Rcpp::traits::input_parameter<const bool>::type adapt_step(adapt_stepSEXP);
  rcpp_result_gen = Rcpp::wrap(solve_final_size_iterative(
      contact_matrix, demography_vector, susceptibility, iterations, tolerance,
      step_rate, adapt_step));
  return rcpp_result_gen;
  END_RCPP
}
// solve_final_size_newton
Eigen::ArrayXd solve_final_size_newton(const Eigen::MatrixXd& contact_matrix,
                                       const Eigen::VectorXd& demography_vector,
                                       const Eigen::VectorXd& susceptibility,
                                       const int iterations,
                                       const double tolerance);
RcppExport SEXP _finalsize_solve_final_size_newton(SEXP contact_matrixSEXP,
                                                   SEXP demography_vectorSEXP,
                                                   SEXP susceptibilitySEXP,
                                                   SEXP iterationsSEXP,
                                                   SEXP toleranceSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Eigen::MatrixXd&>::type contact_matrix(
      contact_matrixSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type demography_vector(
      demography_vectorSEXP);
  Rcpp::traits::input_parameter<const Eigen::VectorXd&>::type susceptibility(
      susceptibilitySEXP);
  Rcpp::traits::input_parameter<const int>::type iterations(iterationsSEXP);
  Rcpp::traits::input_parameter<const double>::type tolerance(toleranceSEXP);
  rcpp_result_gen = Rcpp::wrap(
      solve_final_size_newton(contact_matrix, demography_vector, susceptibility,
                              iterations, tolerance));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_finalsize_epi_spread", (DL_FUNC)&_finalsize_epi_spread, 4},
    {"_finalsize_solve_final_size_iterative",
     (DL_FUNC)&_finalsize_solve_final_size_iterative, 7},
    {"_finalsize_solve_final_size_newton",
     (DL_FUNC)&_finalsize_solve_final_size_newton, 5},
    {NULL, NULL, 0}};

RcppExport void R_init_finalsize(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
