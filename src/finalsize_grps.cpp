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
//' @title Calculate final epidemic size with risk groups using Eigen backend
//'
//' @description This function calculates final epidemic size using SIR model
//' for a heterogeneously mixing population, using an Eigen backend
//'
//' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
//' number of contacts in group $i$ reported by participants in group j
//' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
//' of total population in group $i$ (model will normalise if needed)
//' @param p_susceptibility A matrix giving the probability that an individual
//' in demography group $i$ is in risk (or susceptibility) group $j$.
//' Each row represents the overall distribution of individuals in demographic
//' group $i$ across risk groups, and each row *must sum to 1.0*.
//' @param susceptibility A matrix giving the susceptibility of individuals in
//' demographic group $i$ and risk group $j$.
//' @param solver Which solver to use. Options are "iterative" or "newton", for
//' the iterative solver, or the Newton solver, respectively. Special conditions
//' apply when using the Newton solver.
//' @param iterations Number of solver iterations. Defaults to 1,000.
//' @param tolerance Solver error tolerance.
//' @param step_rate The solver step rate for the iterative solver. Defaults to
//' 1.9 as a value found to work well.
//' @param adapt_step Boolean, whether the iterative solver step rate should be
//' changed based on the solver error. Defaults to TRUE.
//'
//' @keywords epidemic model
//' @export
// [[Rcpp::export]]
Eigen::ArrayXd final_size_grps_cpp(const Eigen::MatrixXd &contact_matrix,
                                   const Eigen::VectorXd &demography_vector,
                                   const Eigen::MatrixXd &p_susceptibility,
                                   const Eigen::MatrixXd &susceptibility,
                                   const Rcpp::String solver,
                                   const int iterations = 1000,
                                   const double tolerance = 1e-6,
                                   double step_rate = 1.9,
                                   const bool adapt_step = true) {
  if (contact_matrix.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: contact matrix must have as many rows as demography groups\n");
  }
  if (p_susceptibility.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: p_susceptibility must have as many rows as demography "
        "groups\n");
  }
  if (susceptibility.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: susceptibility must have as many rows as demography groups\n");
  }
  if (p_susceptibility.size() != susceptibility.size()) {
    Rcpp::stop(
        "Error: p_susceptibility and susceptibility must be matrices of the "
        "same dims\n");
  }
  // check that p_susceptibility rowwise sums are approx 1.0 - ideally sum to 1
  for (size_t i = 0; i < p_susceptibility.rows(); i++) {
    if (std::abs(p_susceptibility.row(i).sum() - 1.0) > 1e-6) {
      Rcpp::stop("Error: p_susceptibility matrix rows must sum to 1.0");
    }
  }
  // check that solver options are correct
  if ((solver != Rcpp::String("iterative")) &&
      (solver != Rcpp::String("newton"))) {
    Rcpp::stop("Error: solver must be one of 'iterative' or 'newton'");
  }

  // prepare epidemiological spread data
  epi_spread_data s = epi_spread(contact_matrix, demography_vector,
                                 p_susceptibility, susceptibility);

  // solve for final sizes (epi_final_size, abbrv to efs_tmp)
  Eigen::ArrayXd efs_tmp;
  if (solver == Rcpp::String("iterative")) {
    efs_tmp = solve_final_size_iterative_cpp(
        s.contact_matrix, s.demography_vector, s.susceptibility, iterations,
        tolerance, step_rate, adapt_step);
  } else {
    efs_tmp =
        solve_final_size_newton_cpp(s.contact_matrix, s.demography_vector,
                                    s.susceptibility, iterations, tolerance);
  }

  // prepare final sizes as matrix with n-demography rows, n-risk-grps cols
  Eigen::Map<Eigen::MatrixXd> efs_tmp_2(
      efs_tmp.data(), demography_vector.size(), p_susceptibility.cols());

  Eigen::MatrixXd psusc = p_susceptibility;
  Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);

  // multiply final sizes from pi_2 with proportions of risk groups -- I think
  Eigen::VectorXd v = efs_tmp_2.array() * lps.array();

  // cast data to n-demography rows, n-risk-grps cols dimensions for return
  Eigen::Map<Eigen::MatrixXd> efs_tmp_3(v.data(), demography_vector.rows(),
                                        p_susceptibility.cols());

  // return row wise sum, one per demo group
  return efs_tmp_3.rowwise().sum();
}
