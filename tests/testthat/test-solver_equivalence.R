# Check for solver equivalence
test_that("Check solver equivalence in final_size_grps_cpp", {
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

  # outcome with iterative solver
  epi_outcome_iterative <- final_size_grps(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "iterative"
  )

  # outcome with newton solver
  epi_outcome_newton <- final_size_grps(
    contact_matrix = r0 * c_matrix,
    demography_vector = d_vector,
    susceptibility = susc,
    p_susceptibility = psusc,
    solver = "newton"
  )

  # check for equivalence in all demography groups
  tolerance <- 1e-6
  invisible(
    Map(
      epi_outcome_iterative, epi_outcome_newton,
      f = function(x, y) {
        expect_equal(x, y, tolerance = tolerance)
      }
    )
  )
})
