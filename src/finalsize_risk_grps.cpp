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
//'
//' @keywords epidemic model
//' @export
// [[Rcpp::export]]
Eigen::ArrayXd final_size_grps_cpp(const Eigen::MatrixXd &contact_matrix,
                                   const Eigen::VectorXd &demography_vector,
                                   const Eigen::MatrixXd &p_susceptibility,
                                   const Eigen::MatrixXd &susceptibility) {
  return solve_final_size_by_susceptibility(contact_matrix, demography_vector,
                                            p_susceptibility, susceptibility);
}
