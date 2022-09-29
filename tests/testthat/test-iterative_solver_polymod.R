# basic test that solver returns numerics within range, for multiple risk groups
test_that("Iterative solver works with polymod data", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
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

  epi_outcome <- solve_final_size_iterative(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]]
  )

  # check that solver returns correct types
  testthat::expect_vector(
    object = epi_outcome,
    ptype = numeric()
  )
  # check that solver returns no nans
  testthat::expect_false(
    any(is.nan(epi_outcome))
  )
  # check that solver returns no nas
  testthat::expect_false(
    any(is.na(epi_outcome))
  )
  # check that solver returns no inf
  testthat::expect_false(
    any(is.infinite(epi_outcome))
  )
  # check that solver returns values within range
  testthat::expect_true(
    all(epi_outcome > 0)
  )
  testthat::expect_true(
    all(epi_outcome < 1)
  )
  # check for size of the vector
  testthat::expect_equal(
    n_demo_grps * n_risk_grps,
    length(epi_outcome)
  )
})
