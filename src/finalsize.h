#ifndef FINALSIZE_H
#define FINALSIZE_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

/// function f1 defined in final_size.R
// f1 <- function(beta2, x) {
//     beta2 %*% (1 - x) + log(x)
//   }
inline Eigen::MatrixXd f1(Eigen::MatrixXd beta2, const Eigen::VectorXd x) {
  Eigen::MatrixXd x_ = (beta2 * ((Eigen::VectorXd::Ones(x.size()) - x))) +
                       ((x.array().log()).matrix());

  return x_;
}

/// function f2 defined in final_size.R
// f2 <- function(beta2, x, size) {
//     -beta2 + diag(size) / x
//   }
inline Eigen::MatrixXd f2(Eigen::MatrixXd beta2, const Eigen::VectorXd x,
                          const size_t &size) {
  // make diagonal matrix of dims [size, size]
  // size should be the number of age groups
  Eigen::VectorXd v = Eigen::VectorXd::Ones(size);
  Eigen::MatrixXd diag_size;
  diag_size.resize(size, size);
  diag_size.fill(0.0);
  diag_size.diagonal() = v;

  // divide diagonal matrix of 1s by x
  // x is the current proportion infected per age group
  for (size_t i = 0; i < diag_size.rows(); i++)
    diag_size.row(i) = diag_size.row(i) / x[i];

  return -beta2 + diag_size;
}

/// function to normalise the demography vector
inline Eigen::VectorXd normalise_demography (const Eigen::VectorXd &demography_vector) {
  return demography_vector / (demography_vector.sum());
}

/// function to get largest real eigenvalue
inline double get_max_real_eigenvalue (const Eigen::MatrixXd &a_matrix) {
  Eigen::EigenSolver<Eigen::MatrixXd> es(a_matrix, false);
  Eigen::MatrixXd eig_vals = es.eigenvalues().real();
  return eig_vals.maxCoeff();
}

/// A function for epidemic spread with susceptibility groups
// taken from Edwin van Leeuwen at https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
inline Rcpp::List epi_spread(const Eigen::MatrixXd &contact_matrix,
                      const Eigen::VectorXd &demography_vector,
                      const Eigen::MatrixXd &p_susceptibility, // risk groups
                      const Eigen::MatrixXd &susceptibility // susc of risk grp?
) {

  // check dimensions

  // count number of risk groups
  auto n_susc_groups = p_susceptibility.cols();
  Eigen::MatrixXd p_susceptibility_ = Eigen::MatrixXd::Ones(p_susceptibility.size(), 1);
  
  // make single column matrix from prop_suscep data,
  // prop_suscep is the prob(suscep) per demography group
  Eigen::MatrixXd psusc = p_susceptibility;
  Eigen::Map<Eigen::MatrixXd> lps(psusc.data(), psusc.size(), 1);

  // replicate demography vector by the number of risk groups
  // and multiply it by the prop_suscep values
  Eigen::ArrayXd demography_vector_ = demography_vector.replicate(n_susc_groups, 1).array() * lps.array();

  Eigen::MatrixXd contact_matrix_ = contact_matrix.replicate(n_susc_groups, n_susc_groups);

  // unroll the risk level matrix
  Eigen::MatrixXd susc = susceptibility;
  Eigen::Map<Eigen::MatrixXd> rm(susc.data(), susc.size(), 1);
  
  return Rcpp::List::create(
    Rcpp::Named("n_risk_groups") = n_susc_groups,
    Rcpp::Named("lps_") = lps,
    Rcpp::Named("contact_matrix_") = contact_matrix_,
    Rcpp::Named("rm_") = rm
  );
}

#endif
