#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' C++ backend to calculate final epidemic size.
//'
// [[Rcpp::export]]
Eigen::ArrayXd final_size_cpp (const double &r0,
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::VectorXd &prop_suscep
) {
    // scale demography vector
    Eigen::VectorXd pp0 = demography_vector / (demography_vector.sum());

    // largest real eigenvalue of the contact matrix
    Eigen::EigenSolver<Eigen::MatrixXd> es(contact_matrix, false);
    Eigen::MatrixXd eig_vals = es.eigenvalues().real();
    double eig_val_max = eig_vals.maxCoeff();
    // scale the next generation matrix for max eigenvalue = r0
    Eigen::MatrixXd mm0 = r0 * (contact_matrix / eig_val_max);
    
    // define transmission matrix A = mm_{ij} * pp_{j} / pp_{i}
    Eigen::Map<Eigen::MatrixXd> pp1(pp0.data(), pp0.size(), 1);
    Eigen::ArrayXd beta1 = (mm0 * prop_suscep).array() / pp1.array();

    return beta1;
}   
