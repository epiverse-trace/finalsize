#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/// function for Newton solver
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
inline Eigen::ArrayXd solve_final_size_newton(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &susceptibility, const int iterations = 10000,
    const double tolerance = 1e-6) {
  // count number of demography groups
  size_t nDim = demography_vector.size();

  Eigen::ArrayXi zeros;  // previously in the settings struct
  zeros.resize(nDim);
  zeros.fill(0);

  Eigen::ArrayXd epi_final_size;  // prev in settings struct
  if (epi_final_size.size() != nDim) {
    epi_final_size.resize(nDim);
    epi_final_size.fill(0.5);
  }

  Eigen::MatrixXd contact_matrix_ = contact_matrix;
  for (size_t i = 0; i < contact_matrix.rows(); ++i) {
    // Check if value should be 0 for (limited) performance increase
    if (demography_vector(i) == 0 || susceptibility(i) == 0 ||
        contact_matrix.row(i).sum() == 0) {
      zeros[i] = 1;
      epi_final_size[i] = 0;
    }
    for (size_t j = 0; j < contact_matrix.cols(); ++j) {
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
  // a function f1 which corresponds to the helper fun f1 in helper_funs.h
  auto f1 = [&contact_matrix_](const Eigen::VectorXd &x,
                               Eigen::VectorXd &&cache) {
    cache =
        contact_matrix_ * (1 - x.array()).matrix() + x.array().log().matrix();
    return std::move(cache);
  };

  Eigen::MatrixXd cache_m = contact_matrix_;
  // a function f2 which corresponds to the helper fun f2 in helper_funs.h
  auto f2 = [&contact_matrix_](const Eigen::VectorXd &x,
                               Eigen::MatrixXd &&cache) {
    cache = (1.0 / x.array()).matrix().asDiagonal();
    cache = -contact_matrix_ + std::move(cache);
    return std::move(cache);
  };

  // a function dx_f that wraps f1, f2, and performs a matrix solve
  auto dx_f = [&f1, &f2](const Eigen::VectorXd &x, Eigen::VectorXd &&cache,
                         Eigen::MatrixXd &&cache_m) {
    cache_m = f2(x, std::move(cache_m));
    cache = -f1(x, std::move(cache));
    cache = cache_m.partialPivLu().solve(std::move(cache));
    return std::move(cache);
  };

  // iterate over n-iterations or until the solver tolerance is met
  Eigen::VectorXd x(nDim);
  x.fill(1e-6);
  double error = 0.0;
  for (auto i = 0; i < iterations; ++i) {
    cache_v = dx_f(x, std::move(cache_v), std::move(cache_m)).array();

    error = cache_v.array().abs().sum();
    x += std::move(cache_v);
    if (error < tolerance) {
      // std::cout << "Iter: " << i << " " << error << std::endl;
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

#endif
