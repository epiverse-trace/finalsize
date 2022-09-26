#include <Rcpp.h>
#include <RcppEigen.h>

#include <iostream>

#include "finalsize.h"

/// Function for final size with susceptibility groups
//' @title Calculate final epidemic size over risk groups with RcppEigen backend
//'
//' @description This function calculates final epidemic size using SIR model
//' for a heterogeneously mixing population, with risk groups
//'
//' @param contact_matrix  Social contact matrix. Entry $mm_{ij}$ gives average
//' number of contacts in group $i$ reported by participants in group j
//' @param demography_vector  Demography vector. Entry $pp_{i}$ gives proportion
//' of total population in group $i$ (model will normalise if needed)
//' @param p_susceptibility WIP.
//' @param susceptibility WIP.
//' @param iterations WIP.
//' @param adapt_step WIP
//' @param tolerance WIP.
//'
//' @keywords epidemic model
//' @export
// [[Rcpp::export]]
Eigen::ArrayXd final_size_grps_cpp(const Eigen::MatrixXd &contact_matrix,
                                   const Eigen::VectorXd &demography_vector,
                                   const Eigen::MatrixXd &p_susceptibility,
                                   const Eigen::MatrixXd &susceptibility,
                                   const int iterations = 1000,
                                   const bool adapt_step = true,
                                   const double tolerance = 1e-6) {
  if (contact_matrix.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: contact matrix must have as many rows as demography groups\n");
  }
  if (p_susceptibility.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: p_susceptibility must have as many rows as demography "
        "groups\n");
  }
  if (susceptibility.rows() != demography_vector.size()) {
    Rcpp::stop(
        "Error: susceptibility must have as many rows as demography groups\n");
  }
  if (p_susceptibility.size() != susceptibility.size()) {
    Rcpp::stop(
        "Error: p_susceptibility and susceptibility must be matrices of the "
        "same dims\n");
  }
  // check that p_susceptibility rowwise sums are approx 1.0 - ideally sum to 1
  for (size_t i = 0; i < p_susceptibility.rows(); i++) {
    if (std::abs(p_susceptibility.row(i).sum() - 1.0) > 1e-6) {
      Rcpp::stop("Error: p_susceptibility matrix rows must sum to 1.0");
    }
  }

  return solve_final_size_by_susceptibility(contact_matrix, demography_vector,
                                            p_susceptibility, susceptibility)
      .rowwise()
      .sum();
}
