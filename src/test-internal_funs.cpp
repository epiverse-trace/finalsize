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

// check that function f1 returns output of expected dimensions
context("Fn f1 return is expected length") {
  test_that("Function f1 returns a 1 row matrix") {
    // set up test inputs
    int n_groups = 3;
    Eigen::MatrixXd beta2;
    beta2.resize(n_groups, n_groups);
    beta2.fill(0.2);  // fill with an arbitrary number

    Eigen::Vector3d x(0.2, 0.2, 0.2);  // same length as n_groups
    CATCH_REQUIRE(f1(beta2, x).size() == n_groups);
    CATCH_REQUIRE(f1(beta2, x).rows() == n_groups);
    CATCH_REQUIRE(f1(beta2, x).cols() == 1);
  }
}

// check that function f2 returns output of expected dimensions
context("Fn f2 return is expected size") {
  test_that("Function f2 returns a 1 row matrix") {
    // set up test inputs
    int n_groups = 3;
    Eigen::MatrixXd beta2;
    beta2.resize(n_groups, n_groups);
    beta2.fill(0.2);  // fill with an arbitrary number

    Eigen::Vector3d x(0.2, 0.2, 0.2);  // same length as n_groups

    CATCH_CHECK(f2(beta2, x, n_groups).size() == beta2.size());
    CATCH_CHECK(f2(beta2, x, n_groups).rows() == n_groups);
    CATCH_CHECK(f2(beta2, x, n_groups).cols() == n_groups);
  }
}

/// check that normalise demography sums to 1.0 and has correct size
context("Fn normalise demography works") {
  Eigen::VectorXd demography(3);
  demography << 1.0, 2.0, 3.0;

  Eigen::VectorXd demography_normalised = normalise_demography(demography);

  test_that("Fn normalise demography sums to 1.0 and has correct size") {
    CATCH_CHECK(demography_normalised.size() == demography.size());
    CATCH_CHECK(demography_normalised.sum() == Approx(1.0).epsilon(1e-6));

    // all elements are less than 1.0
    for (size_t i = 0; i < demography_normalised.size(); i++) {
      CATCH_CHECK(demography_normalised[i] < 1.0);
    }
  }
}

/// check that function to get max real eigenvalue works
context("Fn get_max_real_eigenvalue works") {
  Eigen::MatrixXd a_matrix(2, 2);
  a_matrix << 1.0, 0.0, 0.0, 2.0;

  double max_eigv = get_max_real_eigenvalue(a_matrix);

  test_that("Fn max_real_eigenvalues returns correct eigenvalue") {
    CATCH_CHECK(max_eigv == Approx(2.0).epsilon(1e-6));
  }
}

/// check that scaling of the next-gen matrix is correct
context("Fn scale_nextgen_matrix works") {
  Eigen::MatrixXd a_matrix(2, 2);
  a_matrix << 1.0, 0.0, 0.0, 2.0;

  const double r0 = 1.3;

  Eigen::MatrixXd next_gen_matrix = scale_nextgen_matrix(r0, a_matrix);

  test_that("Fn scale_nextgen_matrix returns correctly scaled matrix") {
    CATCH_CHECK(next_gen_matrix.maxCoeff() == Approx(1.3).epsilon(1e-6));
  }
}
