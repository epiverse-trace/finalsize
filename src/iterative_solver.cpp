// Copyright (c) 2022 finalsize authors. See the LICENCE.md file for more.
#include <Rcpp.h>
#include <RcppEigen.h>

// Enable C++14 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

/// function for iterative solver
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
//' @title Iterative solver for final size.
//'
//' @param contact_matrix Social contact matrix. Entry \eqn{mm_{ij}} gives
//' average number of contacts in group \eqn{i} reported by participants in
//' group \eqn{j}.
//' @param demography_vector Demography vector. Entry \eqn{pp_{i}} gives
//' proportion of total population in group \eqn{i}
//' (model will normalise if needed).
//' @param susceptibility A matrix giving the susceptibility of individuals in
//' demographic group \eqn{i} and risk group \eqn{j}.
//' @param iterations Number of solver iterations. Defaults to 10,000.
//' @param tolerance Solver error tolerance. Solving for final size ends when
//' the error drops below this tolerance. Defaults to set `1e-6`.
//' Larger tolerance values are likely to lead to inaccurate final size
//' estimates.
//' @param step_rate The solver step rate. Defaults to 1.9 as a value found to
//' work well.
//' @param adapt_step Boolean, whether the solver step rate should be changed
//' based on the solver error. Defaults to TRUE.
//' @keywords internal
//'
//' @return A vector of final sizes, of the size (N demography groups *
//' N risk groups).
// [[Rcpp::export(name = ".solve_iterative")]]
Eigen::ArrayXd solve_final_size_iterative(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::VectorXd &susceptibility, const int iterations = 10000,
    const double tolerance = 1e-6, double step_rate = 1.9,
    const bool adapt_step = true) {
  // count number of demography groups
  int nDim = demography_vector.size();

  Eigen::ArrayXi zeros(nDim);
  zeros.fill(0);

  Eigen::ArrayXd epi_final_size(nDim);  // prev in settings struct
  epi_final_size.fill(0.5);

  Eigen::MatrixXd contact_matrix_ = contact_matrix;
  for (int i = 0; i < contact_matrix.rows(); ++i) {
    // Check if value should be 0 for (limited) performance increase
    if (demography_vector(i) == 0 || susceptibility(i) == 0 ||
        contact_matrix.row(i).sum() == 0) {
      zeros[i] = 1;
      epi_final_size[i] = 0;
    }
    for (int j = 0; j < contact_matrix.cols(); ++j) {
      if (zeros[j]) {
        contact_matrix_(i, j) = 0;
      } else {
        // Scale contacts appropriately
        contact_matrix_(i, j) = contact_matrix(i, j) * demography_vector(j);
      }
    }
  }

  Eigen::VectorXd epi_final_size_return(nDim);
  // define functions to minimise error in final size estimate
  auto f = [&contact_matrix_, &susceptibility, &zeros](
               const Eigen::VectorXd &epi_final_size,
               Eigen::VectorXd &epi_final_size_return) {
    Eigen::VectorXd s = contact_matrix_ * (-epi_final_size);
    for (int i = 0; i < contact_matrix_.rows(); ++i) {
      if (zeros[i] == 1) {
        epi_final_size_return[i] = 0;
        continue;
      }
      epi_final_size_return[i] = 1;

      epi_final_size_return(i) -= exp(susceptibility(i) * s(i));
    }
  };

  double current_error = step_rate * nDim;
  double error = NAN;
  const double step_change = 1.1;

  for (auto i = 0; i < iterations; ++i) {
    f(epi_final_size, epi_final_size_return);

    Eigen::ArrayXd dpi = epi_final_size - epi_final_size_return.array();
    error = dpi.abs().sum();
    if (error < tolerance) {
      epi_final_size -= dpi;
      break;
    }

    double change = (current_error - error);
    if (change > 0) {
      epi_final_size -= step_rate * dpi;
      if (adapt_step) {
        step_rate *= step_change;
      }
    } else {
      epi_final_size -= dpi;
      if (adapt_step) step_rate = std::max(0.9 * step_rate, 1.0);
    }
    current_error = error;
  }
  if (current_error / tolerance > 100.0) {
    Rcpp::warning(
        "Solver error > 100x solver tolerance, try increasing iterations");
  }

  // Adjust numerical errors;
  for (auto i = 0; i < epi_final_size.size(); ++i) {
    if (zeros[i] ||
        ((epi_final_size(i) < 0) && (epi_final_size(i) > -tolerance)))
      epi_final_size(i) = 0;
  }
  return epi_final_size;
}
