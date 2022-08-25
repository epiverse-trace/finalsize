#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

namespace finalsize_functions {

  inline Eigen::ArrayXd solve_final_size_newton(
      const Eigen::MatrixXd &contact_matrix, const Eigen::VectorXd &demography,
      const Eigen::VectorXd &susceptibility,
      const bool adapt_step = true,
      const double tolerance = 1e-6) {
    // TODO: accept susceptibility vector only and no p_susceptibility
    // TODO: do we still need settings? --- settings components now passed as args
    size_t nDim = demography.size();

    // TODO: allow this to be cached?
    Eigen::ArrayXi zeros; // previously in the settings struct
    zeros.resize(nDim);
    zeros.fill(0);

    Eigen::ArrayXd pi; // prev in settings struct
    Eigen::ArrayXd pi_return(nDim);
    if (pi.size() != nDim) {
      pi.resize(nDim);
      pi.fill(0.5);
    }

    Eigen::MatrixXd contact_matrixM = contact_matrix;
    for (size_t i = 0; i < contact_matrix.rows(); ++i) {
      // Check if value should be 0 for (limited) performance increase
      if (demography(i) == 0 || susceptibility(i) == 0 ||
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
          contact_matrixM(i, j) = susceptibility(i) * contact_matrix(i, j) * demography(j);
        }
      }
    }

    Eigen::VectorXd cache_v = pi;
    auto f1 = [&contact_matrixM](const Eigen::VectorXd &x, Eigen::VectorXd &&cache) {
      cache = contact_matrixM * (1 - x.array()).matrix() + x.array().log().matrix();
      return std::move(cache);
    };

    Eigen::MatrixXd cache_m = contact_matrixM;
    auto f2 = [&contact_matrixM](const Eigen::VectorXd &x, Eigen::MatrixXd &&cache) {
      cache = (1.0 / x.array()).matrix().asDiagonal();
      cache = -contact_matrixM + std::move(cache);
      return std::move(cache);
    };

    auto dx_f = [&f1, &f2](const Eigen::VectorXd &x, Eigen::VectorXd &&cache,
                          Eigen::MatrixXd &&cache_m) {
      cache_m = f2(x, std::move(cache_m));
      cache = -f1(x, std::move(cache));
      cache = cache_m.partialPivLu().solve(std::move(cache));
      return std::move(cache);
    };

    Eigen::VectorXd x = (1 - pi);
    for (auto i = 0; i < 1000; ++i) {
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

  // struct spreaded_t {
  //   Eigen::MatrixXd p_susc;
  //   Eigen::MatrixXd susc;
  //   Eigen::VectorXd demo_v;
  //   Eigen::MatrixXd m;
  // };

  // spreaded_t epi_spread(const Eigen::MatrixXd &contact_matrix,
  //                       const Eigen::VectorXd &demography,
  //                       const Eigen::MatrixXd &p_susceptibility,
  //                       const Eigen::MatrixXd &susceptibility) {
  //   spreaded_t s;
  //   auto ngroups = p_susceptibility.cols();
  //   s.p_susc = Eigen::MatrixXd::Ones(p_susceptibility.size(), 1);

  //   Eigen::MatrixXd psusc = p_susceptibility;
  //   Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);

  //   s.demo_v = demography.replicate(ngroups, 1).array() * lps.array();

  //   s.m = contact_matrix.replicate(ngroups, ngroups);

  //   Eigen::MatrixXd susc = susceptibility;
  //   Eigen::Map<Eigen::MatrixXd> rm(susc.data(), susc.size(), 1);
  //   s.susc = rm;
  //   return s;
  // }

  // inline Eigen::ArrayXd solve_final_size_by_susceptibility(const Eigen::MatrixXd &contact_matrix,
  //                                         const Eigen::VectorXd &demography,
  //                                         const Eigen::MatrixXd &p_susceptibility,
  //                                         const Eigen::MatrixXd &susceptibility,
  //                                         const bool adapt_step = true,
  //                                         const double tolerance = 1e-6) {
  //   auto s = epi_spread(contact_matrix, demography, p_susceptibility, susceptibility);

  //   Eigen::ArrayXd pi_tmp = solve_final_size_newton(s.m, s.demo_v, s.susc, 
  //     adapt_step, tolerance);
    
  //   Eigen::Map<Eigen::MatrixXd> pi_tmp_2(pi_tmp.data(), demography.rows(),
  //                                 p_susceptibility.cols());
  //   return pi_tmp_2;
  // }
}; // namespace finalsize_functions

//' Final size of an epidemic outbreak.
//' 
//' Final size calculation using a Cpp-based Newton solver.
//' 
// [[Rcpp::export]]
Eigen::ArrayXd solve_final_size_internal(const Eigen::MatrixXd &contact_matrix,
                                                 const Eigen::VectorXd &demography,
                                                 const Eigen::MatrixXd &p_susceptibility,
                                                 const Eigen::MatrixXd &susceptibility,
                                                 const bool adapt_step = true,
                                                 const double tolerance = 1e-6) {
  Eigen::ArrayXd pi_tmp;

  return finalsize_functions::solve_final_size_newton(contact_matrix, demography, susceptibility, adapt_step, tolerance);
  
  // finalsize_functions::solve_final_size_by_susceptibility(contact_matrix, demography, p_susceptibility,
  //                                    susceptibility, adapt_step, tolerance);
  // Eigen::MatrixXd psusc = p_susceptibility;
  // Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);
  // Eigen::VectorXd v = pi_tmp.array() * lps.array();
  // Eigen::Map<Eigen::MatrixXd> pi(v.data(), demography.rows(),
  //                                p_susceptibility.cols());
  // return pi.rowwise().sum();
}