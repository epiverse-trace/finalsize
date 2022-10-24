# Prepare common data
contact_matrix <- matrix(1.0 / 200.0, 2, 2)
demography_vector <- rep(100.0, 2)
psusc <- matrix(1, nrow = 2, ncol = 1)
susc <- psusc

# basic test that solver returns numerics within range
test_that("Final size with iterative solver", {
  # prepare some data for the solver
  r0 <- 1.3

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
    any(is.na(epi_outcome$p_infected))
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
  expect_equal(
    length(demography_vector) * ncol(psusc),
    length(epi_outcome$p_infected)
  )
})

# basic test that solver returns correct answer
test_that("Iterative solver returns correct answer", {
  # prepare some data for the solver
  r0 <- 1.3

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc
  )

  epi_outcome_known <- 1 - exp(-r0 * epi_outcome$p_infected)

  tolerance <- 1e-5
  expect_true(
    all(abs(epi_outcome$p_infected - epi_outcome_known) < tolerance)
  )
})

# basic test that higher r0 results in higher final size
test_that("Iterative solver gives higher final size for larger r0", {
  # prepare some data for the solver
  r0_low <- 1.3
  r0_high <- 3.3

  epi_outcome_low <- final_size(
    r0 = r0_low,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc
  )

  epi_outcome_high <- final_size(
    r0 = r0_high,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc
  )

  expect_true(
    all(epi_outcome_high$p_infected > epi_outcome_low$p_infected)
  )
})
