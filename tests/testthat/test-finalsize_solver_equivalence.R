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

# Check final_size works with Newton solver
test_that("Check finalsize by groups works for Polymod, newton solver", {
  r0 <- 2.0

  n_demo_grps <- length(demography_vector)
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
    tolerance = 1e-6
  )

  # prepare outcome with iterative solver
  epi_outcome_iterative <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative",
    control = control
  )

  # prepare outcome with newton solver
  epi_outcome_newton <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "newton",
    control = control
  )

  # check for equivalence
  expect_equal(
    epi_outcome_iterative$p_infected,
    epi_outcome_newton$p_infected,
    tolerance = 1e-5
  )
})
