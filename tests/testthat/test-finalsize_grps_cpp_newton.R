# Check final_size_grps_cpp works with Newton solver
test_that("Check finalsize by groups works for Polymod, newton solver", {
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

  # Scale contact matrix to demography
  c_matrix <- apply(
    c_matrix, 1, function(r) r / d_vector
  )

  n_demo_grps <- length(d_vector)
  n_risk_grps <- 4

  # prepare p_susceptibility and susceptibility
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
  susc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "newton"
  )

  expect_type(
    epi_outcome, "double"
  )
  # check that values are not NaN
  expect_true(
    all(!is.nan(epi_outcome))
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
    n_demo_grps,
    length(epi_outcome)
  )
})

# Check sensitivity to susceptibility using newton solver
test_that("Check that more susceptible demo-grps have higher final size", {
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

  # Scale contact matrix to demography
  c_matrix <- apply(
    c_matrix, 1, function(r) r / d_vector
  )

  n_demo_grps <- length(d_vector)
  n_risk_grps <- 4

  # prepare p_susceptibility and susceptibility
  # susceptibility of demography group i = 1 is 10x i = 2, 3
  psusc <- matrix(
    data = 1, nrow = n_demo_grps, ncol = n_risk_grps
  )
  psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
  susc <- rbind(
    rep(1, n_risk_grps),
    rep(0.1, n_risk_grps),
    rep(0.1, n_risk_grps)
  )

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "newton"
  )

  expect_vector(
    epi_outcome,
    ptype = numeric()
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
    n_demo_grps,
    length(epi_outcome)
  )
  # check that first group has higher final size
  expect_gt(
    epi_outcome[1], epi_outcome[n_demo_grps]
  )
})

# check for correct final size calculation in complex data case
# using newton solver
test_that("Check final size calculation is correct in complex case", {
  # make a contact matrix
  contact_matrix <- c(
    5.329620e-08, 1.321156e-08, 1.832293e-08, 7.743492e-09, 5.888440e-09,
    2.267918e-09, 1.321156e-08, 4.662496e-08, 1.574182e-08, 1.510582e-08,
    7.943038e-09, 3.324235e-09, 1.832293e-08, 1.574182e-08, 2.331416e-08,
    1.586565e-08, 1.146566e-08, 5.993247e-09, 7.743492e-09, 1.510582e-08,
    1.586565e-08, 2.038011e-08, 1.221124e-08, 9.049331e-09, 5.888440e-09,
    7.943038e-09, 1.146566e-08, 1.221124e-08, 1.545822e-08, 8.106812e-09,
    2.267918e-09, 3.324235e-09, 5.993247e-09, 9.049331e-09, 8.106812e-09,
    1.572736e-08
  ) |> matrix(6, 6)

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

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = r0 * contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susc,
    susceptibility = susc,
    solver = "newton"
  )

  expect_vector(
    epi_outcome,
    ptype = numeric()
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
    length(demography_vector),
    length(epi_outcome)
  )

  # test that final size differs by susceptibility group
  expect_lt(
    epi_outcome[5], epi_outcome[1]
  )
  ratio <- sum(epi_outcome * demography_vector) / sum(demography_vector)
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 0.45)
})
