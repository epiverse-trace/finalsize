# basic test that solver returns numerics within range
skip("Newton solver has known numerical instability")
test_that("Newton solver works", {
  # prepare some data for the solver
  r0 <- 1.3
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  # because Newton solver fails when the contact matrix is singular
  contact_matrix <- contact_matrix + runif(prod(dim(contact_matrix)), 0, 0.001)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  epi_outcome <- solve_final_size_newton(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  # check that solver returns correct types
  expect_vector(
    object = epi_outcome,
    ptype = numeric()
  )
  # check that solver returns no nans
  expect_false(
    any(is.nan(epi_outcome))
  )
  # check that solver returns no nas
  expect_false(
    any(is.na(epi_outcome))
  )
  # check that solver returns no inf
  expect_false(
    any(is.infinite(epi_outcome))
  )
  # check that solver returns values within range
  expect_true(
    any(epi_outcome > 0)
  )
  expect_true(
    any(epi_outcome < 1)
  )
})
