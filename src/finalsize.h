#ifndef FINALSIZE_H
#define FINALSIZE_H

#include <Rcpp.h>
#include <RcppEigen.h>

#include "helper_funs.h"

// [[Rcpp::plugins(cpp11)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

/// a struct to hold intermediate outputs
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
struct epi_spread_data {
  Eigen::MatrixXd contact_matrix;
  Eigen::VectorXd demography_vector;
  Eigen::MatrixXd p_susceptibility;
  Eigen::MatrixXd susceptibility;
};

/// function to scale next generation matrix
inline Eigen::MatrixXd scale_nextgen_matrix(
    const double &r0, const Eigen::MatrixXd &contact_matrix) {
  const double max_real_eigv = get_max_real_eigenvalue(contact_matrix);
  return r0 * (contact_matrix / max_real_eigv);
}

/// a struct to hold intermediate outputs
struct epi_spread_data {
  Eigen::MatrixXd contact_matrix;
  Eigen::VectorXd demography_vector;
  Eigen::MatrixXd p_susceptibility;
  Eigen::MatrixXd susceptibility;
};

/// A function for epidemic spread with susceptibility groups
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
inline epi_spread_data epi_spread(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &p_susceptibility,  // risk groups
    const Eigen::MatrixXd &susceptibility     // susc of risk grp?
) {
  // count number of risk groups
  const int n_susc_groups = p_susceptibility.cols();
  Eigen::MatrixXd p_susceptibility_ =
      Eigen::MatrixXd::Ones(p_susceptibility.size(), 1);

  // make single column matrix from prop_suscep data,
  // prop_suscep is the prob(suscep) per demography group
  Eigen::MatrixXd psusc = p_susceptibility;
  Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);

  // replicate demography vector by the number of risk groups
  // and multiply it by the prop_suscep values
  Eigen::VectorXd demography_vector_ =
      demography_vector.replicate(n_susc_groups, 1).array() * lps.array();

  Eigen::MatrixXd contact_matrix_ =
      contact_matrix.replicate(n_susc_groups, n_susc_groups);

  // unroll the risk level matrix
  Eigen::MatrixXd susc = susceptibility;
  Eigen::Map<Eigen::MatrixXd> rm(susc.data(), susc.size(), 1);

  epi_spread_data tmp_data;
  tmp_data.contact_matrix = contact_matrix_;
  tmp_data.demography_vector = demography_vector_;
  tmp_data.p_susceptibility = p_susceptibility_;
  tmp_data.susceptibility = rm;

  return tmp_data;
}

/// function for Newton solver
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
inline Eigen::ArrayXd solve_final_size_newton(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &p_susceptibility,
    const Eigen::MatrixXd &susceptibility, 
    const int iterations = 1000,
    const bool adapt_step = true, const double tolerance = 1e-6) {
  
  // count number of demography groups
  size_t nDim = demography_vector.size();

  Eigen::ArrayXi zeros;  // previously in the settings struct
  zeros.resize(nDim);
  zeros.fill(0);

  Eigen::ArrayXd pi;  // prev in settings struct
  if (pi.size() != nDim) {
    pi.resize(nDim);
    pi.fill(0.5);
  }

  Eigen::MatrixXd contact_matrixM = contact_matrix;
  for (size_t i = 0; i < contact_matrix.rows(); ++i) {
    // Check if value should be 0 for (limited) performance increase
    if (demography_vector(i) == 0 || susceptibility(i) == 0 ||
        contact_matrix.row(i).sum() == 0) {
      zeros[i] = 1;
      pi[i] = 0;
    }
    for (size_t j = 0; j < contact_matrix.cols(); ++j) {
      if (zeros[j] == 1) {
        contact_matrixM(i, j) = 0;
      } else {
        // Scale contacts appropriately
        // Could add transmissibility (j)?
        contact_matrixM(i, j) =
            susceptibility(i) * contact_matrix(i, j) * demography_vector(j);
      }
    }
  }

  Eigen::VectorXd cache_v = pi;
  // a function f1 which corresponds to the helper fun f1 in helper_funs.h
  auto f1 = [&contact_matrixM](const Eigen::VectorXd &x,
                               Eigen::VectorXd &&cache) {
    cache =
        contact_matrixM * (1 - x.array()).matrix() + x.array().log().matrix();
    return std::move(cache);
  };

  Eigen::MatrixXd cache_m = contact_matrixM;
  // a function f2 which corresponds to the helper fun f2 in helper_funs.h
  auto f2 = [&contact_matrixM](const Eigen::VectorXd &x,
                               Eigen::MatrixXd &&cache) {
    cache = (1.0 / x.array()).matrix().asDiagonal();
    cache = -contact_matrixM + std::move(cache);
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

  // iterate over n-iterations (1000) or until the solver tolerance is met
  Eigen::VectorXd x = (1 - pi);
  for (auto i = 0; i < iterations; ++i) {
    cache_v = dx_f(x, std::move(cache_v), std::move(cache_m)).array();

    double error = cache_v.array().abs().sum();
    x += std::move(cache_v);
    if (error < tolerance) {
      // std::cout << "Iter: " << i << " " << error << std::endl;
      break;
    }
  }

  pi = 1 - x.array();

  return pi;
}

/// function to solve final size by susceptibility
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
inline Eigen::MatrixXd solve_final_size_by_susceptibility(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &p_susceptibility,
    const Eigen::MatrixXd &susceptibility, 
    const int iterations = 1000,
    const bool adapt_step = true, const double tolerance = 1e-6) {

  // prepare epidemiological spread data
  epi_spread_data s = epi_spread(contact_matrix, demography_vector,
                                 p_susceptibility, susceptibility);

  // solve for final sizes
  Eigen::ArrayXd pi =
      solve_final_size_newton(s.contact_matrix, s.demography_vector,
                              s.p_susceptibility, s.susceptibility);

  // prepare final sizes as matrix with n-demography rows, n-risk-grps cols
  Eigen::Map<Eigen::MatrixXd> pi_2(pi.data(), demography_vector.size(),
                                   p_susceptibility.cols());

  Eigen::MatrixXd psusc = p_susceptibility;
  Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);

  // multiply final sizes from pi_2 with proportions of risk groups -- I think
  Eigen::VectorXd v = pi_2.array() * lps.array();

  // cast data to n-demography rows, n-risk-grps cols dimensions for return
  Eigen::Map<Eigen::MatrixXd> pi_3(v.data(), demography_vector.rows(),
                                   p_susceptibility.cols());

  return pi_3;
}

#endif
