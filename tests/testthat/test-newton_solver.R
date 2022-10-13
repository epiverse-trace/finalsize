# basic test that solver returns numerics within range
test_that("Final size with newton solver", {
  # prepare some data for the solver
  r0 <- 1.3
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  epi_outcome <- final_size(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
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
    all(epi_outcome > 0)
  )
  expect_true(
    all(epi_outcome < 1)
  )
  # check for size of the vector
  expect_equal(
    length(demography_vector) * ncol(psusc),
    length(epi_outcome)
  )
})

# basic test that newton solver returns correct answer
test_that("Newton solver returns correct answer", {
  # prepare some data for the solver
  r0 <- 1.3
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  epi_outcome <- final_size(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
  )

  epi_outcome_known <- 1 - exp(-r0 * epi_outcome)

  tolerance <- 1e-5
  expect_true(
    all(abs(epi_outcome - epi_outcome_known) < tolerance)
  )
})

# basic test that higher r0 results in higher final size
test_that("Newton solver gives higher final size for large r0", {
  # prepare some data for the solver
  r0_low <- 1.3
  r0_high <- 3.3
  contact_matrix <- matrix(1 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  epi_outcome_low <- final_size(
    contact_matrix = r0_low * contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
  )

  epi_outcome_high <- final_size(
    contact_matrix = r0_high * contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
  )

  expect_true(
    all(epi_outcome_high > epi_outcome_low)
  )
})
