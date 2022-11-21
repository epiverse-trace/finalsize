# Calculates the upper limit of final size given the r0
# The upper limit is given by a well mixed population
upper_limit <- function(r0) {
  f <- function(par) {
    abs(1 - exp(-r0 * par[1]) - par[1])
  }
  opt <- optim(
    par = 0.5, fn = f,
    lower = 0, upper = 1,
    method = "Brent"
  )
  opt
}

# Prepare common data
contact_matrix <- matrix(1.0 / 200.0, 2, 2)
demography_vector <- rep(100.0, 2)
psusc <- matrix(1, nrow = 2, ncol = 1)
susc <- psusc

# basic test that solver works for various r0 near 2.0 and above
test_that("Iterative solver works with r0 = 2", {
  # prepare some data for the solver
  r0 <- 2

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
    any(epi_outcome$p_infected > 0)
  )
  expect_true(
    any(epi_outcome$p_infected < 1)
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector) * ncol(psusc)
  )
})

# check solver works for r0 = 4
test_that("Iterative solver works with r0 = 4", {
  # prepare some data for the solver
  r0 <- 4

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
  # check for correct answer
  tolerance <- 1e-5
  expected_outcome <- rep(
    upper_limit(r0)$par,
    length(demography_vector)
  )
  expect_equal(
    epi_outcome$p_infected, expected_outcome,
    tolerance = tolerance
  )

  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector) * ncol(psusc)
  )
})

# check solver works for r0 = 12
test_that("Iterative solver works with r0 = 12", {
  # prepare some data for the solver
  r0 <- 12

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    control = list(
      tolerance = 1e-3
    )
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
  # check for correct answer
  tolerance <- 1e-5
  expected_outcome <- rep(
    upper_limit(r0)$par,
    length(demography_vector)
  )
  expect_equal(
    epi_outcome$p_infected, expected_outcome,
    tolerance = tolerance
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector) * ncol(psusc)
  )
})

# check solver works for locked adaptive step size
test_that("Iterative solver works with r0 = 4, locked step size", {
  # prepare some data for the solver
  r0 <- 4
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  psusc <- matrix(1, nrow = 2, ncol = 1)
  susc <- psusc

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    control = list(
      adapt_step = FALSE
    )
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
  # check for correct answer
  expect_equal
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector) * ncol(psusc)
  )
})
