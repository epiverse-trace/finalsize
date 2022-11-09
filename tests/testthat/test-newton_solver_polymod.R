# Prepare common elements for testing
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- contact_data$matrix
demography_vector <- contact_data$demography$population

# scale by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- contact_matrix / demography_vector

# basic test that solver returns numerics within range, for multiple risk groups
test_that("Newton solver works with polymod data", {
  r0 <- 1.3

  n_demo_grps <- length(demography_vector)
  n_risk_grps <- 3

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- psusc / rowSums(psusc)
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
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

# Newton solver works with Polymod data with r0 = 2
test_that("Newton solver works with polymod data", {
  r0 <- 2.0

  n_demo_grps <- length(demography_vector)
  n_risk_grps <- 3

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- psusc / rowSums(psusc)
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
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
