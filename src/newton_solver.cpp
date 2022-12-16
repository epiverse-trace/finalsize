// Copyright (c) 2022 finalsize authors. See the LICENCE.md file for more.
#include <Rcpp.h>
#include <RcppEigen.h>

// Enable C++14 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

/// function for Newton solver
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
//' @title Newton solver for final size.
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
//' @keywords internal
//'
//' @return A two dimensional array of final sizes per age-risk group.
// [[Rcpp::export(name = ".solve_newton")]]
Eigen::ArrayXd solve_final_size_newton(const Eigen::MatrixXd &contact_matrix,
                                       const Eigen::VectorXd &demography_vector,
                                       const Eigen::VectorXd &susceptibility,
                                       const int iterations = 10000,
                                       const double tolerance = 1e-6) {
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
      if (zeros[j] == 1) {
        contact_matrix_(i, j) = 0;
      } else {
        // Scale contacts appropriately
        // Could add transmissibility (j)?
        contact_matrix_(i, j) =
            susceptibility(i) * contact_matrix(i, j) * demography_vector(j);
      }
    }
  }

  Eigen::VectorXd cache_v = epi_final_size;
  // a function f1 that multiplies the contact matrix by the final size guess
  // + the log of the guess
  auto f1 = [&contact_matrix_](const Eigen::VectorXd &x,
                               Eigen::VectorXd &cache) {
    cache = -(contact_matrix_ * (1 - x.array()).matrix() +
              x.array().log().matrix());
  };

  Eigen::MatrixXd cache_m = contact_matrix_;
  // a function f2 which adds the negative of the contact matrix
  // to a diagonal matrix of the current final size guess
  auto f2 = [&contact_matrix_](const Eigen::VectorXd &x,
                               Eigen::MatrixXd &cache) {
    cache = (1.0 / x.array()).matrix().asDiagonal();
    cache = -contact_matrix_ + cache;
  };

  // a function dx_f that wraps f1, f2, and performs a matrix solve
  auto dx_f = [&f1, &f2](const Eigen::VectorXd &x, Eigen::VectorXd &cache,
                         Eigen::MatrixXd &cache_m) {
    f2(x, cache_m);
    f1(x, cache);
    cache = cache_m.partialPivLu().solve(cache).array();
  };

  // iterate over n-iterations or until the solver tolerance is met
  Eigen::VectorXd x(nDim);
  x.fill(1e-6);
  double error = 0.0;
  for (auto i = 0; i < iterations; ++i) {
    dx_f(x, cache_v, cache_m);

    error = cache_v.array().abs().sum();
    x += cache_v;
    if (error < tolerance) {
      break;
    }
  }
  if (error / tolerance > 100.0) {
    Rcpp::warning(
        "Solver error > 100x solver tolerance, try increasing iterations");
  }

  epi_final_size = 1 - x.array();

  return epi_final_size;
}
