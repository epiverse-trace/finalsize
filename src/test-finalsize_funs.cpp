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

/// check that fun epi_spread returns an epi_spread_data object w/ correct attrs
context("Fn epi_spread works for a single risk group") {
  // make some function arguments
  Eigen::MatrixXd contact_matrix(2, 2);
  contact_matrix << 1.0, 0.5, 0.5, 1.0;
  Eigen::VectorXd demography_vector(2);
  demography_vector << 1.0, 1.0;

  // make p_susceptibility with one (1) risk group
  Eigen::MatrixXd p_susceptibility(2, 1);  // n rows per demo grp, 1 col risk
  p_susceptibility.col(0) << 1.0, 1.0;

  Eigen::MatrixXd susceptibility = p_susceptibility;

  epi_spread_data tmp_data = epi_spread(contact_matrix, demography_vector,
                                        p_susceptibility, susceptibility);

  test_that("Fn epi_spread returns contaxt matrix as is") {
    CATCH_CHECK(tmp_data.contact_matrix.size() == contact_matrix.size());
    CATCH_CHECK(tmp_data.contact_matrix.rows() == contact_matrix.rows());
    CATCH_CHECK(tmp_data.contact_matrix.cols() == contact_matrix.cols());
  }

  test_that("Fn epi_spread returns demography vector as is") {
    CATCH_CHECK(tmp_data.demography_vector.size() == demography_vector.size());
  }
}

/// check that fun epi_spread works for multiple risk groups
context("Fn epi_spread works for a multiple risk groups") {
  // make some function arguments
  Eigen::MatrixXd contact_matrix(2, 2);
  contact_matrix << 1.0, 0.5, 0.5, 1.0;
  Eigen::VectorXd demography_vector(2);
  demography_vector << 1.0, 1.0;

  // make p_susceptibility with one (1) risk group
  Eigen::MatrixXd p_susceptibility(2, 2);  // n rows per demo grp, 2 col risk
  p_susceptibility << 0.25, 0.75, 0.25, 0.75;

  Eigen::MatrixXd susceptibility = p_susceptibility;
  susceptibility.fill(1.0);  // all groups have same susceptibility?

  epi_spread_data tmp_data = epi_spread(contact_matrix, demography_vector,
                                        p_susceptibility, susceptibility);

  test_that("Fn epi_spread returns replicated contaxt matrix") {
    // CATCH_CHECK(tmp_data.contact_matrix.size() == contact_matrix.size());
    CATCH_CHECK(tmp_data.contact_matrix.rows() ==
                (contact_matrix.rows() * p_susceptibility.cols()));
    CATCH_CHECK(tmp_data.contact_matrix.cols() ==
                (contact_matrix.cols() * p_susceptibility.cols()));
  }

  test_that("Fn epi_spread returns replicated demography vector") {
    CATCH_CHECK(tmp_data.demography_vector.size() ==
                (demography_vector.size() * p_susceptibility.cols()));
  }

  test_that("Fn epi_spread stores p_susceptibility size but not dims, vals") {
    CATCH_CHECK(tmp_data.p_susceptibility.size() == p_susceptibility.size());
    CATCH_CHECK(tmp_data.p_susceptibility.rows() == p_susceptibility.size());
    CATCH_CHECK(tmp_data.p_susceptibility.cols() == 1);
  }
}

/// check Newton solver
context("Newton solver works") {
  // make some function arguments
  Eigen::MatrixXd contact_matrix(2, 2);
  contact_matrix << 1.0, 0.5, 0.5, 1.0;
  Eigen::VectorXd demography_vector(2);
  demography_vector << 1.0, 1.0;

  // make p_susceptibility with one (1) risk group
  Eigen::MatrixXd p_susceptibility(2, 2);  // n rows per demo grp, 2 col risk
  p_susceptibility << 0.25, 0.75, 0.25, 0.75;

  Eigen::MatrixXd susceptibility = p_susceptibility;
  susceptibility.fill(1.0);  // all groups have same susceptibility?

  epi_spread_data tmp_data = epi_spread(contact_matrix, demography_vector,
                                        p_susceptibility, susceptibility);

  Eigen::ArrayXd final_size_temp = solve_final_size_newton(
      contact_matrix, demography_vector, p_susceptibility, susceptibility);

  test_that("Fn solve_final_size_newton returns same size as demography") {
    CATCH_CHECK(final_size_temp.size() == demography_vector.size());
  }
  test_that("Fn solve_final_size_newton values LTE 1.0") {
    for (size_t i = 0; i < final_size_temp.size(); i++) {
      CATCH_CHECK(final_size_temp[i] <= 1.0);
    }
  }
}

/// check that solving by susceptibility groups works
context("Solving by susceptibility groups works") {
  // make some function arguments
  Eigen::MatrixXd contact_matrix(2, 2);
  contact_matrix << 1.0, 0.5, 0.5, 1.0;
  Eigen::VectorXd demography_vector(2);
  demography_vector << 1.0, 1.0;

  // make p_susceptibility with one (1) risk group
  Eigen::MatrixXd p_susceptibility(2, 2);  // n rows per demo grp, 2 col risk
  p_susceptibility << 0.25, 0.75, 0.25, 0.75;

  Eigen::MatrixXd susceptibility = p_susceptibility;
  susceptibility.fill(1.0);

  Eigen::MatrixXd final_size = solve_final_size_by_susceptibility(
      contact_matrix, demography_vector, p_susceptibility, susceptibility, 1001,
      true, 1.1e-6);

  test_that("Solving by susceptibility returns same rows as demography") {
    CATCH_CHECK(final_size.rows() == demography_vector.size());
  }
  test_that("Solving by susceptibility size is demography X risk groups") {
    CATCH_CHECK(final_size.size() ==
                (demography_vector.size() * p_susceptibility.cols()));
  }
}
