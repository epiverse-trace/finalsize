#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

inline Eigen::MatrixXd f1 (Eigen::MatrixXd &beta2, const Eigen::VectorXd &x) {
    
    Eigen::VectorXd x_ = ((Eigen::VectorXd::Ones(x.size()) - x)) + 
        ((x.array().log()).matrix());
    
    for (size_t i = 0; i < beta2.rows(); i++)
        beta2.row(i) = beta2.row(i) * x_[i];
    
    return beta2;
}
//' C++ backend to calculate final epidemic size.
//'
// [[Rcpp::export]]
Rcpp::List final_size_cpp (const double &r0,
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const double &prop_suscep
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
    Eigen::MatrixXd beta1 = (mm0 * prop_suscep).array();
    // rowwise division by age group proportion
    // this is vectorised in final_size.R
    // with demography vector recycled over columns; i.e., age 1 * row 1 etc.
    for (size_t i = 0; i < beta1.rows(); i++)
        beta1.row(i) = beta1.row(i) / pp0[i];

    Eigen::MatrixXd beta2 = beta1.transpose();
    // rowwise multiplication with corresponding age group proportion
    for (size_t i = 0; i < beta2.rows(); i++)
        beta2.row(i) = beta2.row(i) * pp0[i];
    
    beta2 = (beta2.array().transpose()).matrix(); // could be more concise?

    // Newton solver for final size equation A(1-x) = -log(x)
    // get number of age groups
    const size_t v_size = pp0.size();

    // define functions f1 and f2
    
    return Rcpp::List::create(
        pp0,
        mm0,
        prop_suscep,
        beta1,
        beta2,
        v_size,
        f1 (beta2, pp0)
    );
}
