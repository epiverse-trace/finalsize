#pragma once

#include "finalsize.h"

#include <Rcpp.h>
#include <RcppEigen.h>

#include <iostream>

// [[Rcpp::plugins(cpp11)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' @title Calculate final epidemic size with RcppEigen backend
//'
//' @description This function calculates final epidemic size using SIR model
//' for a heterogeneously mixing population.
//'
//' @param r0  Basic reproduction number.
//' @param contact_matrix  Social contact matrix. Entry $mm_{ij}$ gives average
//' number of contacts in group $i$ reported by participants in group j
//' @param demography_vector  Demography vector. Entry $pp_{i}$ gives proportion
//' of total population in group $i$ (model will normalise if needed)
//' @param prop_initial_infected Proportion of all age groups that is initially
//' infected. May be a single number, or a vector of proportions infected.
//' If a vector, must be the same length as the demography vector, otherwise the
//' vector will be recycled. Default value is 0.001 for all groups.
//' @param prop_suscep  Proportion of each group susceptible. May be a single
//' numeric value or a numeric vector of the same length as the demography
//' vector.
//' @param iterations Number of solver iterations
//'
//' @keywords epidemic model
//' @export
// [[Rcpp::export]]
Eigen::VectorXd final_size_cpp(const double &r0,
                               const Eigen::MatrixXd &contact_matrix,
                               const Eigen::VectorXd &demography_vector,
                               const Eigen::VectorXd &prop_initial_infected,
                               Eigen::VectorXd prop_suscep,
                               const int iterations = 30) {
  if (prop_suscep.size() == 1) {
    const double prop_suscep_ = prop_suscep[0];
    prop_suscep.resize(demography_vector.size(), 1);
    prop_suscep.fill(prop_suscep_);
    Rcpp::Rcout << "prop_suscep = " << prop_suscep << "\n";
  } else if (prop_suscep.size() != demography_vector.size()) {
    Rcpp::stop("Error: prop_suscep must be same size as demography vector");
  }

  // scale demography vector
  Eigen::VectorXd pp0 = demography_vector / (demography_vector.sum());

  // largest real eigenvalue of the contact matrix
  Eigen::EigenSolver<Eigen::MatrixXd> es(contact_matrix, false);
  Eigen::MatrixXd eig_vals = es.eigenvalues().real();
  double eig_val_max = eig_vals.maxCoeff();
  // scale the next generation matrix for max eigenvalue = r0
  Eigen::MatrixXd mm0 = r0 * (contact_matrix / eig_val_max);

  // define transmission matrix A = mm_{ij} * pp_{j} / pp_{i}
  Eigen::MatrixXd beta1 = mm0.array();
  // rowwise division by age group proportion
  // this is vectorised in final_size.R
  // with demography vector recycled over columns; i.e., age 1 * row 1 etc.
  for (size_t i = 0; i < beta1.rows(); i++)
    beta1.row(i) = beta1.row(i) * prop_suscep[i] / pp0[i];

  Eigen::MatrixXd beta2 = beta1.transpose();
  // rowwise multiplication with corresponding age group proportion
  for (size_t i = 0; i < beta2.rows(); i++)
    beta2.row(i) = beta2.row(i) * pp0[i];

  beta2 = (beta2.array().transpose()).matrix();  // could be more concise?

  // Newton solver for final size equation A(1-x) = -log(x)
  // get number of age groups
  const size_t v_size = pp0.size();

  // prepare timesteps as iterations matrix
  // fill with -1.0
  Eigen::MatrixXd iterate_output;
  iterate_output.resize(static_cast<size_t>(iterations), v_size);
  iterate_output.fill(-1.0);

  // prepare initial infections vector
  Eigen::VectorXd x0;
  if (prop_initial_infected.size() == 1) {
    x0 = prop_initial_infected[0] * pp0;
  } else {
    if (prop_initial_infected.size() != pp0.size()) {
      Rcpp::stop(
          "Error: prop_initial_infection must be same size as demography "
          "vector\n");
    }
    Rcpp::Rcout << "multiple prop_initial_infected\n";

    x0 = prop_initial_infected.array() * pp0.array();
  }

  // iterate over timesteps
  iterate_output.row(0) = x0;  // initial proportion infected

  // holding matrix (vec) for dx, the change in infection proportions
  Eigen::MatrixXd f1_m, f2_m, dx;
  // update the size of the current outbreak
  // take the 1st (0) column of dx, as all cols per row are same
  // iterate_output.row(i) = iterate_output.row(i - 1) + dx.col(0).array();
  // starting at second row
  for (size_t i = 1; i < static_cast<size_t>(iterations); i++) {
    f2_m = f2(beta2, iterate_output.row(i - 1), v_size);
    f1_m = -f1(beta2, iterate_output.row(i - 1));
    dx = f2_m.partialPivLu().solve(f1_m);
    // crude iteration-and-column wise addition
    for (size_t j = 0; j < iterate_output.cols(); j++) {
      iterate_output(i, j) = iterate_output(i - 1, j) + dx.col(0)[j];
    }
  }

  // return 1 - (vector of final proportions of each age group infected)
  return (1.0 - (iterate_output.row(iterate_output.rows() - 1)).array());
}
