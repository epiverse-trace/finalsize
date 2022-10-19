#ifndef EPI_SPREAD_H
#define EPI_SPREAD_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/// a struct to hold intermediate outputs
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
struct epi_spread_data {
  Eigen::MatrixXd contact_matrix;
  Eigen::VectorXd demography_vector;
  Eigen::VectorXd susceptibility;
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

  // make single column matrix from prop_suscep data,
  // prop_suscep is the prob(suscep) per demography group
  Eigen::MatrixXd psusc = p_susceptibility;
  const Eigen::VectorXd lps(
      Eigen::Map<Eigen::VectorXd>(psusc.data(), psusc.size()));

  // replicate demography vector by the number of risk groups
  // and multiply it by the prop_suscep values
  const Eigen::VectorXd demography_vector_ =
      demography_vector.replicate(n_susc_groups, 1).array() * lps.array();

  const Eigen::MatrixXd contact_matrix_ =
      contact_matrix.replicate(n_susc_groups, n_susc_groups);

  // unroll the risk level matrix
  Eigen::MatrixXd susc = susceptibility;
  const Eigen::VectorXd rm(
      Eigen::Map<Eigen::VectorXd>(susc.data(), susc.size()));

  epi_spread_data tmp_data;
  tmp_data.contact_matrix = contact_matrix_;
  tmp_data.demography_vector = demography_vector_;
  tmp_data.susceptibility = rm;

  return tmp_data;
}

#endif
