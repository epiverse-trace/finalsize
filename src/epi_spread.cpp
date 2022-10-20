#ifndef EPI_SPREAD_H
#define EPI_SPREAD_H

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/// A function for epidemic spread with susceptibility groups
// taken from Edwin van Leeuwen at
// https://gitlab.com/epidemics-r/code_snippets/feature/newton_solver/include/finalsize.hpp
//' @title Prepare population data for solvers.
//'
//' @description An internal function that prepares the contact matrix,
//' demography data, and susceptibility matrix for solvers.
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
//'
//' @return A list object with named elements: "contact_matrix",
//' "demography_vector", "p_susceptibility_", and "susceptibility".
//' The contact matrix is replicated row and column wise for each risk group
//' and the demography vector is replicated for each risk group.
// [[Rcpp::export(name = ".epi_spread")]]
Rcpp::List epi_spread(
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

  return Rcpp::List::create(
      Rcpp::Named("contact_matrix") = contact_matrix_,
      Rcpp::Named("demography_vector") = demography_vector_,
      Rcpp::Named("susceptibility") = rm);
}

#endif
