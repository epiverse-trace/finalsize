/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <Rcpp.h>
#include <RcppEigen.h>
#include <testthat.h>

#include "finalsize.h"

/// some example data
// These were generated using:
// cm_lst <- contactsr::load_polymod_data(age_limits=c(15,30,45,60,75))
// contactsr::symmetric_matrix(cm_lst)
// taken from
// https://gitlab.com/epidemics-r/code_snippets/-/blob/feature/newton_solver/tests/catch_final_size.cpp
auto example_contact_matrix() {
  Eigen::Matrix<double, 6, 6> cm;
  cm << 5.329620e-08, 1.321156e-08, 1.832293e-08, 7.743492e-09, 5.888440e-09,
      2.267918e-09, 1.321156e-08, 4.662496e-08, 1.574182e-08, 1.510582e-08,
      7.943038e-09, 3.324235e-09, 1.832293e-08, 1.574182e-08, 2.331416e-08,
      1.586565e-08, 1.146566e-08, 5.993247e-09, 7.743492e-09, 1.510582e-08,
      1.586565e-08, 2.038011e-08, 1.221124e-08, 9.049331e-09, 5.888440e-09,
      7.943038e-09, 1.146566e-08, 1.221124e-08, 1.545822e-08, 8.106812e-09,
      2.267918e-09, 3.324235e-09, 5.993247e-09, 9.049331e-09, 8.106812e-09,
      1.572736e-08;
  return cm;
}

auto example_demography() {
  Eigen::Array<double, 6, 1> demo;
  demo << 10831795, 11612456, 13511496, 11499398, 8167102, 4587765;
  return demo;
}

/// check that solving by risk groups returns the correct answer
context("Correct solution to final size") {
  const double r0 = 1.3;
  Eigen::MatrixXd contact_matrix(2, 2);
  contact_matrix.fill(r0 / 200.0);
  Eigen::ArrayXd demography_vector(2);
  demography_vector.fill(100.0);

  Eigen::ArrayXd p_susceptibility(2);
  p_susceptibility.fill(1);

  Eigen::ArrayXd susceptibility(2);
  susceptibility.fill(1);

  Eigen::MatrixXd pi = solve_final_size_by_susceptibility(
      contact_matrix, demography_vector, p_susceptibility, susceptibility);

  Eigen::ArrayXd pi_2 = pi.rowwise().sum();
  test_that("Correct solutions to final size simple case") {
    CATCH_CHECK(pi_2(0) > 0.0);
    CATCH_CHECK(pi_2(0) == pi_2(1));
    CATCH_CHECK(Approx(pi_2(0)) == 1 - exp(-r0 * pi_2(0)));
  }
}

/// check that final size by groups cpp for polymod matrix
context("Correct solutions for Polymod matrix") {
  double r0 = 1.3;
  auto lambda = example_contact_matrix();
  auto demography = example_demography();

  Eigen::ArrayXd psusc(lambda.rows());
  psusc.fill(1);

  Eigen::ArrayXd susc(lambda.rows());
  susc.fill(1);

  Eigen::VectorXd pi =
      solve_final_size_by_susceptibility(r0 * lambda, demography, psusc, susc);

  test_that("Correct solutions for Polymod matrix") {
    for (size_t i = 0; i < pi.size(); i++) {
      CATCH_CHECK(pi(i) > 0.0);
    }

    CATCH_CHECK(pi(5) < pi(0));

    double ratio = (pi.array() * demography).sum() / demography.sum();
    CATCH_CHECK(ratio > 0.3);
    CATCH_CHECK(ratio < 0.45);
  }
}
