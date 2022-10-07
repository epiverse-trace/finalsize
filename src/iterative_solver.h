#ifndef ITERATIVE_SOLVER_H
#define ITERATIVE_SOLVER_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/// function for iterative solver
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
Eigen::ArrayXd solve_final_size_iterative_cpp(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &p_susceptibility,
    const Eigen::MatrixXd &susceptibility, const int iterations = 1000,
    const double tolerance = 1e-6, double step_rate = 1.9,
    const bool adapt_step = true) {
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

  Eigen::VectorXd epi_final_size_return(nDim);
  // define functions to minimise error in final size estimate
  auto f = [&contact_matrix_, &p_susceptibility, &susceptibility, &zeros](
               const Eigen::VectorXd &epi_final_size,
               Eigen::VectorXd &&epi_final_size_return) {
    Eigen::VectorXd s = contact_matrix_ * (-epi_final_size);
    for (size_t i = 0; i < contact_matrix_.rows(); ++i) {
      if (zeros[i] == 1) {
        epi_final_size_return[i] = 0;
        continue;
      }
      epi_final_size_return[i] = 1;

      for (size_t k = 0; k < p_susceptibility.cols(); ++k) {
        epi_final_size_return(i) -=
            p_susceptibility(i, k) * exp(susceptibility(i, k) * s(i));
      }
    }
    return std::move(epi_final_size_return);
  };

  double current_error = step_rate * nDim;

  for (auto i = 0; i < 10000; ++i) {
    epi_final_size_return = f(epi_final_size, std::move(epi_final_size_return));

    Eigen::ArrayXd dpi = epi_final_size - epi_final_size_return.array();
    double error = dpi.abs().sum();
    if (error < tolerance) {
      epi_final_size -= dpi;
      break;
    }

    double change = (current_error - error);
    if (change > 0) {
      epi_final_size -= step_rate * dpi;
      if (adapt_step) step_rate *= 1.1;
    } else {
      epi_final_size -= dpi;
      if (adapt_step) step_rate = std::max(0.9 * step_rate, 1.0);
    }
    current_error = error;
  }

  // Adjust numerical errors;
  for (auto i = 0; i < epi_final_size.size(); ++i) {
    if (zeros[i])
      // if (epi_final_size(i) < 0 && epi_final_size(i) > -tolerance)
      epi_final_size(i) = 0;
  }
  return epi_final_size;
}

#endif
