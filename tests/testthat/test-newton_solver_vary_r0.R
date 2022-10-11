# Calculates the upper limit of final size given the r0
# The upper limit is given by a well mixed population
upper_limit <- function(r0) {
  f <- function(par) {
    abs(1 - exp(-r0 * par[1]) - par[1])
  }
  opt <- optim(
    par = c(0.5), fn = f,
    lower = c(0), upper = c(1),
    method = "Brent"
  )
  opt
}

# basic test that Newton solver works for various r0 near 2.0 and above
test_that("Newton solver works with r0 = 2", {
  # prepare some data for the solver
  r0 <- 2
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
    any(epi_outcome > 0)
  )
  expect_true(
    any(epi_outcome < 1)
  )
  # check for size of the vector
  expect_equal(
    length(demography_vector) * ncol(psusc),
    length(epi_outcome)
  )
})

# check Newton solver works for r0 = 4
test_that("Newton solver works with r0 = 4", {
  # prepare some data for the solver
  r0 <- 4
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
  # check for correct answer
  tolerance <- 1e-5
  lapply(epi_outcome, \(x) {
    expect_equal(
      x, upper_limit(r0)$par,
      tolerance = tolerance
    )
  })
  # check for size of the vector
  expect_equal(
    length(demography_vector) * ncol(psusc),
    length(epi_outcome)
  )
})

# check Newton solver works for r0 = 12
test_that("Newton solver works with r0 = 12", {
  # prepare some data for the solver
  r0 <- 12
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
  # check for correct answer
  tolerance <- 1e-5
  lapply(epi_outcome, \(x) {
    expect_equal(
      x, upper_limit(r0)$par,
      tolerance = tolerance
    )
  })
  # check for size of the vector
  expect_equal(
    length(demography_vector) * ncol(psusc),
    length(epi_outcome)
  )
})
