# basic test that solver returns numerics within range, for multiple risk groups
test_that("Newton solver works with polymod data", {
  r0 <- 1.3
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40),
    split = TRUE
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion

  n_demo_grps <- length(d_vector)
  n_risk_grps <- 3

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
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
    n_demo_grps * n_risk_grps,
    length(epi_outcome)
  )
})

# Newton solver works with Polymod data with r0 = 2
test_that("Newton solver works with polymod data", {
  r0 <- 2.0
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40),
    split = TRUE
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion

  n_demo_grps <- length(d_vector)
  n_risk_grps <- 3

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  # needed to get demography-risk combinations
  epi_spread_data <- epi_spread(
    contact_matrix = c_matrix,
    demography_vector = d_vector,
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
    n_demo_grps * n_risk_grps,
    length(epi_outcome)
  )
})
