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
Rcpp::NumericVector final_size_cpp (const double &r0,
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector
) {
    // scale demography vector
    Eigen::VectorXd pp0 = demography_vector / (demography_vector.sum());

    // largest real eigenvalue of the contact matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(contact_matrix, false);
    double eigv_max = es.eigenvalues().maxCoeff();

    return Rcpp::NumericVector::create(eigv_max);
}
