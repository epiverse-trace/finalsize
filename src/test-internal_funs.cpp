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
#include <testthat.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "finalsize.hpp"

// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.
int twoPlusTwo() {
  return 2 + 2;
}

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Sample unit tests") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("two plus two equals four") {
    expect_true(twoPlusTwo() == 4);
  }

}

context("Fn f1 return is expected length") {
  
  test_that("Function f1 returns a 1 row matrix") {
    // set up test inputs
    int n_groups = 3;
    Eigen::MatrixXd beta2;
    beta2.resize(n_groups, n_groups);
    beta2.fill(0.2); // fill with an arbitrary number

    Eigen::Vector3d x(0.2, 0.2, 0.2); // same length as n_groups
    CATCH_REQUIRE(f1(beta2, x).size() == n_groups);
    CATCH_REQUIRE(f1(beta2, x).rows() == n_groups);
    CATCH_REQUIRE(f1(beta2, x).cols() == 1);
  }
}

context("Fn f2 return is expected size") {
  
  test_that("Function f1 returns a 1 row matrix") {
    // set up test inputs
    int n_groups = 3;
    Eigen::MatrixXd beta2;
    beta2.resize(n_groups, n_groups);
    beta2.fill(0.2); // fill with an arbitrary number

    Eigen::Vector3d x(0.2, 0.2, 0.2); // same length as n_groups

    CATCH_REQUIRE(f2(beta2, x, n_groups).size() == beta2.size());
    CATCH_REQUIRE(f2(beta2, x, n_groups).rows() == n_groups);
    CATCH_REQUIRE(f2(beta2, x, n_groups).cols() == n_groups);
  }
}
