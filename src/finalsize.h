#ifndef FINALSIZE_H
#define FINALSIZE_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

/// a struct to hold intermediate outputs
struct epi_spread_data {
  Eigen::MatrixXd contact_matrix;
  Eigen::VectorXd demography_vector;
  Eigen::MatrixXd p_susceptibility;
  Eigen::MatrixXd susceptibility;
};

inline Eigen::MatrixXd scale_nextgen_matrix(
    const double &r0, const Eigen::MatrixXd &contact_matrix) {
  const double max_real_eigv = get_max_real_eigenvalue(contact_matrix);
  return r0 * (contact_matrix / max_real_eigv);
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

  epi_spread_data tmp_data;
  tmp_data.contact_matrix = contact_matrix_;
  tmp_data.demography_vector = demography_vector_;
  tmp_data.p_susceptibility = p_susceptibility_;
  tmp_data.susceptibility = rm;

  return tmp_data;
}

#endif
