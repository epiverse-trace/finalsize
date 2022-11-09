# basic test that solver returns numerics within range, for multiple risk groups
test_that("Iterative solver works with multiple risk groups", {
  r0 <- 1.3
  # prepare some data for the solver
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)

  # p_susceptibility and susceptibility
  n_risk_grps <- 3
  psusc <- matrix(1, nrow = 2, ncol = n_risk_grps)
  psusc <- psusc / rowSums(psusc) # rows sum to 1.0
  susc <- matrix(1, nrow = length(demography_vector), ncol = n_risk_grps)

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc
  )

  # check that solver returns correct types
  expect_s3_class(
    object = epi_outcome,
    "data.frame"
  )
  # check that solver returns no nans
  expect_false(
    any(is.nan(epi_outcome$p_infected))
  )
  # check that solver returns no nas
  expect_false(
    anyNA(epi_outcome$p_infected)
  )
  # check that solver returns no inf
  expect_false(
    any(is.infinite(epi_outcome$p_infected))
  )
  # check that solver returns values within range
  expect_true(
    all(epi_outcome$p_infected > 0)
  )
  expect_true(
    all(epi_outcome$p_infected < 1)
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector) * n_risk_grps
  )
})

# olver returns numerics within range, for multiple risk groups and demo groups
test_that("Iterative solver works with multiple risk and age groups", {
  r0 <- 1.3
  # prepare some data for the solver
  n_demo_grps <- 5
  contact_matrix <- matrix(r0 / 200.0, nrow = n_demo_grps, ncol = n_demo_grps)
  demography_vector <- rep(100.0, n_demo_grps)

  # p_susceptibility and susceptibility
  n_risk_grps <- 3
  psusc <- matrix(1, nrow = n_demo_grps, ncol = n_risk_grps)
  psusc <- psusc / rowSums(psusc) # rows sum to 1.0
  susc <- matrix(1, nrow = n_demo_grps, ncol = n_risk_grps)

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc
  )

  # check that solver returns correct types
  expect_s3_class(
    object = epi_outcome,
    "data.frame"
  )
  # check that solver returns no nans
  expect_false(
    any(is.nan(epi_outcome$p_infected))
  )
  # check that solver returns no nas
  expect_false(
    anyNA(epi_outcome$p_infected)
  )
  # check that solver returns no inf
  expect_false(
    any(is.infinite(epi_outcome$p_infected))
  )
  # check that solver returns values within range
  expect_true(
    all(epi_outcome$p_infected > 0)
  )
  expect_true(
    all(epi_outcome$p_infected < 1)
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    n_demo_grps * n_risk_grps
  )
})
