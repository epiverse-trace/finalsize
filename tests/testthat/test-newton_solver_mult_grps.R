# basic test that Newton solver returns within range, for multiple risk groups
test_that("Newton solver works with multiple risk groups", {
  r0 <- 1.3
  # prepare some data for the solver
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)

  # p_susceptibility and susceptibility
  n_risk_grps <- 3
  psusc <- matrix(1, nrow = 2, ncol = n_risk_grps)
  psusc <- t(apply(psusc, 1, \(x) x / sum(x))) # rows sum to 1.0
  susc <- matrix(1, nrow = length(demography_vector), ncol = n_risk_grps)

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome <- solve_final_size_newton(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
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
    length(demography_vector) * n_risk_grps,
    length(epi_outcome)
  )
})

# Newton solver returns within range, for multiple risk groups and demo groups
test_that("Newton solver works with multiple risk and age groups", {
  r0 <- 1.3
  # prepare some data for the solver
  demo_grps <- 5
  contact_matrix <- matrix(r0 / 200.0, nrow = demo_grps, ncol = demo_grps)
  demography_vector <- rep(100.0, demo_grps)

  # p_susceptibility and susceptibility
  n_risk_grps <- 3
  psusc <- matrix(1, nrow = demo_grps, ncol = n_risk_grps)
  psusc <- t(apply(psusc, 1, \(x) x / sum(x))) # rows sum to 1.0
  susc <- matrix(1, nrow = demo_grps, ncol = n_risk_grps)

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome <- solve_final_size_newton(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
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
    demo_grps * n_risk_grps,
    length(epi_outcome)
  )
})
