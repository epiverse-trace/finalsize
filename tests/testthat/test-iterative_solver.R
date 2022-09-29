# basic test that solver returns numerics within range
test_that("Iterative solver works", {
  # prepare some data for the solver
  r0 <- 1.3
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome <- solve_final_size_iterative(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]]
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

# basic test that solver returns correct answer
test_that("Iterative solver returns correct answer", {
  # prepare some data for the solver
  r0 <- 1.3
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome <- solve_final_size_iterative(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]]
  )

  epi_outcome_known <- 1 - exp(-r0 * epi_outcome)

  tolerance <- 1e-5
  expect_true(
    all(abs(epi_outcome - epi_outcome_known) < tolerance)
  )
})

# basic test that higher r0 results in higher final size
test_that("Iterative solver returns correct answer", {
  # prepare some data for the solver
  r0_low <- 1.3
  r0_high <- 3.3
  contact_matrix <- matrix(1 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome_low <- solve_final_size_iterative(
    contact_matrix = r0_low * epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]]
  )

  epi_outcome_high <- solve_final_size_iterative(
    contact_matrix = r0_high * epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]]
  )

  expect_true(
    all(epi_outcome_high > epi_outcome_low)
  )
})
