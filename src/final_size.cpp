#include <Rcpp.h>
#include <RcppEigen.h>

#include <iostream>

#include "epi_spread.h"
#include "iterative_solver.h"
#include "newton_solver.h"

// [[Rcpp::plugins(cpp11)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

/// Function for final size with susceptibility groups
//' @title Final size of an epidemic
//'
//' @description Calculates the final size of an epidemic outbreak
//' in a population with heterogeneous mixing, and with heterogeneous
//' susceptibility to infection such as that conferred by an immunisation
//' programme.
//'
//' **Note**: This is an internal function that is called by
//' \link[finalsize]{final_size}. This function is somewhat faster, but lacks
//' argument checks found in `final_size`. Use with caution.
//'
//' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives the
//' average number of contacts in group \eqn{i} reported by participants in
//' group \eqn{j}.
//' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
//' proportion of total population in group \eqn{i} (model will normalise
//' if needed).
//' @param p_susceptibility A matrix giving the probability that an individual
//' in demography group \eqn{i} is in risk (or susceptibility) group \eqn{j}.
//' Each row represents the overall distribution of individuals in demographic
//' group \eqn{i} across risk groups, and each row *must sum to 1.0*.
//' @param susceptibility A matrix giving the susceptibility of individuals in
//' demographic group \eqn{i} and risk group \eqn{j}.
//' @param solver Which solver to use. Options are "iterative" or "newton", for
//' the iterative solver, or the Newton solver, respectively. The Newton solver
//' only uses the `iterations` and `tolerance` options.
//' @param iterations The number of iterations over which to solve for the final
//' size, unless the error is below the solver tolerance.
//' @param tolerance The solver tolerance, set to `1e-6` by default; solving for
//' final size ends when the error drops below this tolerance.
//' @param step_rate The solver step rate. Defaults to 1.9 as a value found to
//' work well.
//' @param adapt_step Boolean, whether the solver step rate should be changed
//' based on the solver error. Defaults to TRUE.
//'
//' @keywords epidemic model
// [[Rcpp::export(name = ".final_size")]]
Eigen::ArrayXd final_size_(const Eigen::MatrixXd &contact_matrix,
                           const Eigen::VectorXd &demography_vector,
                           const Eigen::MatrixXd &p_susceptibility,
                           const Eigen::MatrixXd &susceptibility,
                           const Rcpp::String &solver = "iterative",
                           const int &iterations = 10000,
                           const double &tolerance = 1e-6,
                           const double &step_rate = 1.9,
                           const bool &adapt_step = true) {
  // prepare epidemiological spread data
  epi_spread_data s = epi_spread(contact_matrix, demography_vector,
                                 p_susceptibility, susceptibility);

  // pointer to solver function, iterative by default
  Eigen::ArrayXd efs_tmp;
  if (solver == "iterative") {
    efs_tmp = solve_final_size_iterative(s.contact_matrix, s.demography_vector,
                                         s.susceptibility, iterations,
                                         tolerance, step_rate, adapt_step);
  } else if (solver == "newton") {
    efs_tmp = solve_final_size_newton(s.contact_matrix, s.demography_vector,
                                      s.susceptibility, iterations, tolerance);
  } else {
    Rcpp::stop("Error: solver must be one of 'iterative' or 'newton'");
  }

  // multiply final sizes from pi_2 with proportions of risk groups -- I think
  efs_tmp = efs_tmp * s.p_susceptibility.array();

  // cast data to n-demography rows, n-risk-grps cols dimensions for return
  const Eigen::MatrixXd epi_final_size(Eigen::Map<Eigen::MatrixXd>(
      efs_tmp.data(), demography_vector.rows(), p_susceptibility.cols()));

  // return row wise sum, one per demo group
  return epi_final_size.rowwise().sum();
}
