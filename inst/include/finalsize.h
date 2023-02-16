// Copyright 2023 'finalsize' authors. See repository licence in LICENSE.md.
// Written manually to allow header export
// copied from https://github.com/r-pkg-examples/rcpp-shared-cpp-functions

#ifndef finalsize_finalsize_H_
#define finalsize_finalsize_H_

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

// Include the Rcpp Header and RcppEigen
#include <Rcpp.h>
#include <RcppEigen.h>

#include "iterative_solver.h"  // paths must be relative to this file
#include "newton_solver.h"

#endif  // finalsize_finalsize_H_
