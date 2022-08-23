#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>

//' Cpp calculation of final outbreak size.
//'
//' @description Test function to calculate squares.
//'
//' @param contact_matrix Social contact matrix. Entry mm_ij gives average
//' number of contacts in group i reported by participants in group j.
//' @param demography Demography vector. Entry pp_i gives proportion of
//' total population in group i (model will normalise if needed)
//' @param susceptibility  Proportion of each group susceptible. Null assumption is
//' fully susceptible.
//' @param adapt_step Solver step size is adaptive? WIP.
//' @param step_rate Solver step rate. WIP
//' @param cache_pi Whether to cache matrix? WIP.
//' @param tolerance WIP.
//' @keywords epidemic model
//' @return Final size of the epidemic as a vector, with one element per age and susceptibility group.
//' @export
// [[Rcpp::export]] 
Eigen::ArrayXd final_size_newton_cpp(
    const Eigen::MatrixXd &contact_matrix, const Eigen::VectorXd &demography,
    const Eigen::VectorXd &susceptibility,
    const bool &adapt_step = true,
    const double &step_rate = 1.9,
    const bool &cache_pi = true,
    const double &tolerance = 1e-6) {
    // TODO: accept susceptibility vector only and no p_susceptibility
    // TODO: do we still need settings? --- settings removed and passed as arguments.
    size_t nDim = demography.size();

    // TODO: allow this to be cached?
    Eigen::ArrayXi zeros;
    zeros.resize(nDim);
    zeros.fill(0);

    // an unitialised array called pi;
    Eigen::ArrayXd pi;
    Eigen::ArrayXd pi_return(nDim);
    if (pi.size() != nDim || !cache_pi) {
        pi.resize(nDim);
        pi.fill(0.5);
    }

    Eigen::MatrixXd lambdaM = contact_matrix;
    for (size_t i = 0; i < contact_matrix.rows(); ++i) {
    // Check if value should be 0 for (limited) performance increase
        if (demography(i) == 0 || susceptibility(i) == 0 ||
            contact_matrix.row(i).sum() == 0) {
            zeros[i] = 1;
            pi[i] = 0;
        }
        for (size_t j = 0; j < contact_matrix.cols(); ++j) {
            if (zeros[j] == 1) {
                lambdaM(i, j) = 0;
            } else {
                // Scale contacts appropriately
                // Could add transmissibility (j)?
                lambdaM(i, j) = susceptibility(i) * contact_matrix(i, j) * demography(j);
            }
        }
    }

    Eigen::VectorXd cache_v = pi;
    auto f1 = [&lambdaM](const Eigen::VectorXd &x, Eigen::VectorXd &&cache) {
        cache = lambdaM * (1 - x.array()).matrix() + x.array().log().matrix();
        return std::move(cache);
    };

    Eigen::MatrixXd cache_m = lambdaM;
    auto f2 = [&lambdaM](const Eigen::VectorXd &x, Eigen::MatrixXd &&cache) {
        cache = (1.0 / x.array()).matrix().asDiagonal();
        cache = -lambdaM + std::move(cache);
        return std::move(cache);
    };

    auto dx_f = [&f1, &f2](const Eigen::VectorXd &x, Eigen::VectorXd &&cache,
                            Eigen::MatrixXd &&cache_m) {
        cache_m = f2(x, std::move(cache_m));
        cache = -f1(x, std::move(cache));
        cache = cache_m.partialPivLu().solve(std::move(cache));
        return std::move(cache);
    };

    Eigen::VectorXd x = (1 - pi);
    for (auto i = 0; i < 1000; ++i) {
        cache_v = dx_f(x, std::move(cache_v), std::move(cache_m)).array();

        double error = cache_v.array().abs().sum();
        x += std::move(cache_v);
        if (error < tolerance) {
            // std::cout << "Iter: " << i << " " << error << std::endl;
            break;
        }
    }
    pi = 1 - x.array();
    return pi;
}


