# Prepare common testing elements
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

# Check final_size works with iterative solver
test_that("Check finalsize by groups works for Polymod, iterative solver", {
  r0 <- 2.0
  n_demo_grps <- length(demography_vector)
  n_risk_grps <- 4L

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- psusc / rowSums(psusc)
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  # prepare control
  control <- list(
    iterations = 10000,
    tolerance = 1e-6,
    step_rate = 1.9,
    adapt_step = TRUE
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative",
    control = control
  )

  expect_s3_class(
    epi_outcome, "data.frame"
  )
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
    all(epi_outcome$p_infected > 0)
  )
  expect_true(
    all(epi_outcome$p_infected < 1)
  )
  # check for number of demography groups returned
  expect_identical(
    n_demo_grps,
    length(unique(epi_outcome$demo_grp))
  )
  # check for overall length
  expect_identical(
    n_demo_grps * n_risk_grps,
    nrow(epi_outcome)
  )

  # expect error if solver option is incorrect
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = psusc,
      susceptibility = susc,
      solver = "some wrong option",
      control = list()
    ),
    regexp = "(Error)*(should be one of)*(iterative)*(newton)"
  )
})

# Check sensitivity to susceptibility using newton solver
test_that("Check that more susceptible demo-grps have higher final size", {
  r0 <- 1.3

  n_demo_grps <- length(demography_vector)
  n_risk_grps <- 4

  # prepare p_susceptibility and susceptibility
  # susceptibility of demography group i = 1 is 10x i = 2, 3
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- psusc / rowSums(psusc)
  susc <- rbind(
    rep(1, n_risk_grps),
    rep(0.1, n_risk_grps),
    rep(0.1, n_risk_grps)
  )

  # prepare control
  control <- list(
    iterations = 10000,
    tolerance = 1e-6,
    step_rate = 1.9,
    adapt_step = TRUE
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative",
    control = control
  )

  expect_s3_class(
    epi_outcome,
    "data.frame"
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
  # check that first group has higher final size
  expect_gt(
    epi_outcome$p_infected[1], epi_outcome$p_infected[n_demo_grps]
  )
})

# check for correct final size calculation in complex data case
# using newton solver
test_that("Check final size calculation is correct in complex case", {
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
  susc <- p_susc

  # prepare control
  control <- list(
    iterations = 10000,
    tolerance = 1e-6,
    step_rate = 1.9,
    adapt_step = TRUE
  )

  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susc,
    susceptibility = susc,
    solver = "iterative",
    control = control
  )

  expect_s3_class(
    epi_outcome,
    "data.frame"
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
