// Copyright 2023 'finalsize' authors. See repository licence in LICENSE.md.
#include <Rcpp.h>
#include <finalsize.h>
#include <string.h>

//' @title Calculate the final size of an epidemic
//' @description An internal function that interfaces between the R function
//' `final_size()` and functions in the package header.
//' @param parameters A named list of parameters for the final size calculation.
//' See the R function documentation for details and input checking.
//'
// [[Rcpp::export(name=".final_size")]]
Eigen::ArrayXd final_size(const Rcpp::List &parameters) {
  if (strcmp(parameters["solver"], "iterative") == 0) {
    return finalsize::solve_final_size_iterative(
        parameters["contact_matrix"], parameters["demography_vector"],
        parameters["susceptibility"], parameters["iterations"],
        parameters["tolerance"], parameters["step_rate"],
        parameters["adapt_step"]);
  } else {
    return finalsize::solve_final_size_newton(
        parameters["contact_matrix"], parameters["demography_vector"],
        parameters["susceptibility"], parameters["iterations"],
        parameters["tolerance"]);
  }
}
