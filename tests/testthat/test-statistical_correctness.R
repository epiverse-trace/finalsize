#### Test statistical correctness of final_size output ####

# Calculates the upper limit of final size given the r0
# The upper limit is given by a well mixed population
upper_limit <- function(r0) {
  f <- function(par) {
    abs(1 - exp(-r0 * par[1]) - par[1])
  }
  opt <- optim(
    par = 0.5, fn = f,
    lower = 0.0, upper = 1.0,
    method = "Brent"
  )
  opt
}

#### Prepare synthetic data ####
contact_matrix <- matrix(1.0 / 200.0, nrow = 2L, ncol = 2L)
demography_vector <- rep(100.0, 2L)
p_susceptibility <- matrix(1.0, nrow = 2L, ncol = 1L)
susceptibility <- p_susceptibility

#### Test that final_size outputs are correct and valid across R0 ####

# Test final size outputs are valid using the iterative solver
test_that("final_size iterative solver outputs are valid for R0 values", {
  r0_values <- c(1.9, 2.0, 4.0, 12.0)

  # store success final_size p_infected values
  p_infected_list <- list()

  for (r0 in r0_values) {
    epi_outcome <- final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "iterative",
      control = list(
        iterations = 100000, # Increased for tests
        tolerance = 1e-4 # Reduced for R0 = 12.0
      )
    )

    p_infected_list[[as.character(r0)]] <- epi_outcome$p_infected

    # check that values are not NaN
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
      all(epi_outcome$p_infected >= 0.0)
    )
    expect_true(
      all(epi_outcome$p_infected <= 1.0)
    )

    # check for correct answer
    tolerance <- 1e-4
    expected_outcome <- rep(
      upper_limit(r0)$par,
      length(demography_vector)
    )
    expect_equal(
      epi_outcome$p_infected, expected_outcome,
      tolerance = tolerance,
      info = sprintf("Failed for R0 = %f", r0)
    )
  }

  # Check that larger R0 leads to larger final size
  expect_true(
    all(p_infected_list[[1]] < p_infected_list[[4]])
  )
})

#### Test that final_size outputs are correct and valid across R0 ####

# Test final size outputs are valid using the Newton solver
test_that("final_size Newton solver outputs are valid for R0 values", {
  r0_values <- c(1.9, 2.0, 4.0, 12.0)

  # store success final_size p_infected values
  p_infected_list <- list()

  for (r0 in r0_values) {
    epi_outcome <- final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "newton",
      control = list(
        iterations = 100000, # Increased for tests
        tolerance = 1e-5 # Newton solver requires lower tolerance
      )
    )

    p_infected_list[[as.character(r0)]] <- epi_outcome$p_infected

    # check that values are not NaN
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
      all(epi_outcome$p_infected >= 0.0)
    )
    expect_true(
      all(epi_outcome$p_infected <= 1.0)
    )

    # check for correct answer
    tolerance <- 1e-4
    expected_outcome <- rep(
      upper_limit(r0)$par,
      length(demography_vector)
    )
    expect_equal(
      epi_outcome$p_infected, expected_outcome,
      tolerance = tolerance,
      info = sprintf("Failed for R0 = %f", r0)
    )
  }

  # Check that larger R0 leads to larger final size
  expect_true(
    all(p_infected_list[[1]] < p_infected_list[[4]])
  )
})

#### Check for solver equivalence with synthetic data ####
# Test for equivalence across R0 values
test_that("Solvers return equivalent solutions with synthetic data", {
  r0_values <- c(1.9, 2.0, 4.0, 12.0)

  for (r0 in r0_values) {
    epi_outcome_iterative <- final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "iterative",
      control = list(
        tolerance = 1e-4 # Newton solver requires lower tolerance
      )
    )

    epi_outcome_newton <- final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "newton",
      control = list(
        tolerance = 1e-5 # Newton solver requires lower tolerance
      )
    )

    # check for correct answer
    tolerance <- 1e-4
    expect_equal(
      epi_outcome_iterative$p_infected, epi_outcome_newton$p_infected,
      tolerance = tolerance,
      info = sprintf("Failed for R0 = %f", r0)
    )
  }
})

#### Check solver equivalence with POLYMOD data ####

# Prepare common elements for testing; POLYMOD data
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

r0 <- 2.0

n_demo_grps <- length(demography_vector)
n_risk_grps <- 3L

# prepare p_susceptibility and susceptibility
p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)

# multiple susceptibility groups
susceptibility <- matrix(
  data = seq(0.1, 1.0, length.out = n_risk_grps),
  nrow = n_demo_grps, ncol = n_risk_grps,
  byrow = TRUE
)

# assign column names
colnames(susceptibility) <- c("immunised", "part-immunised", "susceptible")

# prepare control
control <- list(
  iterations = 10000,
  tolerance = 1e-6
)

# prepare outcome with iterative solver
epi_outcome_iterative <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "iterative",
  control = control
)

# prepare outcome with newton solver
epi_outcome_newton <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "newton",
  control = control
)

# Check final_size works with Newton solver
test_that("Solvers return equivalent solutions with POLYMOD data", {
  # check for equivalence
  expect_equal(
    epi_outcome_iterative$p_infected,
    epi_outcome_newton$p_infected,
    tolerance = 1e-5
  )
})

### Test that lower susceptibility leads to lower final size ####

test_that("Lower susceptibility leads to lower final size", {
  # all fully susceptibles must have larger final size than immunised
  expect_true(
    all(epi_outcome_newton[epi_outcome_newton$susc_grp ==
      "susceptible", ]$p_infected >
      epi_outcome_newton[epi_outcome_newton$susc_grp ==
        "immunised", ]$p_infected)
  )
})


#### Check for correct final size calculation in complex data case ####

# test taken from EvL
test_that("Newton solver is correct in complex case", {
  # make a contact matrix
  contact_matrix <- matrix(
    data = c(
      5.329620e-08, 1.321156e-08, 1.832293e-08, 7.743492e-09, 5.888440e-09,
      2.267918e-09, 1.321156e-08, 4.662496e-08, 1.574182e-08, 1.510582e-08,
      7.943038e-09, 3.324235e-09, 1.832293e-08, 1.574182e-08, 2.331416e-08,
      1.586565e-08, 1.146566e-08, 5.993247e-09, 7.743492e-09, 1.510582e-08,
      1.586565e-08, 2.038011e-08, 1.221124e-08, 9.049331e-09, 5.888440e-09,
      7.943038e-09, 1.146566e-08, 1.221124e-08, 1.545822e-08, 8.106812e-09,
      2.267918e-09, 3.324235e-09, 5.993247e-09, 9.049331e-09, 8.106812e-09,
      1.572736e-08
    ),
    nrow = 6, ncol = 6
  )

  # make a demography vector
  demography_vector <- c(
    10831795, 11612456, 13511496,
    11499398, 8167102, 4587765
  )

  # get an example r0
  r0 <- 1.3

  # a p_susceptibility matrix
  p_susc <- matrix(1, nrow(contact_matrix), 1)
  susceptibility <- p_susc

  # prepare control
  control <- list(
    iterations = 10000,
    tolerance = 1e-6
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susc,
    susceptibility = susceptibility,
    solver = "newton",
    control = control
  )

  # check that solver returns values within range
  expect_true(
    all(epi_outcome$p_infected >= 0)
  )
  expect_true(
    all(epi_outcome$p_infected <= 1)
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector)
  )

  # test that final size differs by susceptibility group
  expect_lt(
    epi_outcome$p_infected[5], epi_outcome$p_infected[1]
  )
  ratio <- sum(epi_outcome$p_infected * demography_vector) /
    sum(demography_vector)
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 0.45)
})

#### Check for correct final size calculation in complex data case ####

# test taken from EvL
test_that("Iterative solver is correct in complex case", {
  # make a contact matrix
  contact_matrix <- matrix(
    data = c(
      5.329620e-08, 1.321156e-08, 1.832293e-08, 7.743492e-09, 5.888440e-09,
      2.267918e-09, 1.321156e-08, 4.662496e-08, 1.574182e-08, 1.510582e-08,
      7.943038e-09, 3.324235e-09, 1.832293e-08, 1.574182e-08, 2.331416e-08,
      1.586565e-08, 1.146566e-08, 5.993247e-09, 7.743492e-09, 1.510582e-08,
      1.586565e-08, 2.038011e-08, 1.221124e-08, 9.049331e-09, 5.888440e-09,
      7.943038e-09, 1.146566e-08, 1.221124e-08, 1.545822e-08, 8.106812e-09,
      2.267918e-09, 3.324235e-09, 5.993247e-09, 9.049331e-09, 8.106812e-09,
      1.572736e-08
    ),
    nrow = 6, ncol = 6
  )

  # make a demography vector
  demography_vector <- c(
    10831795, 11612456, 13511496,
    11499398, 8167102, 4587765
  )

  # get an example r0
  r0 <- 1.3

  # a p_susceptibility matrix
  p_susc <- matrix(1, nrow(contact_matrix), 1)
  susceptibility <- p_susc

  # prepare control
  control <- list(
    iterations = 10000,
    tolerance = 1e-6
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susc,
    susceptibility = susceptibility,
    solver = "iterative",
    control = control
  )

  # check that solver returns values within range
  expect_true(
    all(epi_outcome$p_infected >= 0)
  )
  expect_true(
    all(epi_outcome$p_infected <= 1)
  )
  # check for size of the vector
  expect_length(
    epi_outcome$p_infected,
    length(demography_vector)
  )

  # test that final size differs by susceptibility group
  expect_lt(
    epi_outcome$p_infected[5], epi_outcome$p_infected[1]
  )
  ratio <- sum(epi_outcome$p_infected * demography_vector) /
    sum(demography_vector)
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 0.45)
})
