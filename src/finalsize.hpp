#pragma once

/// function f1 defined in final_size.R
// f1 <- function(beta2, x) {
//     beta2 %*% (1 - x) + log(x)
//   }
inline Eigen::MatrixXd f1 (Eigen::MatrixXd beta2, const Eigen::VectorXd x) {
    
    Eigen::MatrixXd x_ = (beta2 * ((Eigen::VectorXd::Ones(x.size()) - x))) + 
        ((x.array().log()).matrix());
    
    return x_;
}

/// function f2 defined in final_size.R
// f2 <- function(beta2, x, size) {
//     -beta2 + diag(size) / x
//   }
inline Eigen::MatrixXd f2 (Eigen::MatrixXd beta2, const Eigen::VectorXd x,
    const size_t &size
) {
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
