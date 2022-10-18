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

  # prepare outcome with iterative solver
  epi_outcome_iterative <- final_size_grps_cpp(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative",
    control = control
  )

  # prepare outcome with newton solver
  epi_outcome_newton <- final_size_grps_cpp(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "newton",
    control = control
  )

  # check for equivalence
  expect_equal(epi_outcome_iterative, epi_outcome_newton, tolerance = 1e-5)
})
