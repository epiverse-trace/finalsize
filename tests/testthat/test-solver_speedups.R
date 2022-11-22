#### Tests minor solver options ####

#### Prepare data ####
r0 <- 2.0
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- contact_data$matrix
demography_vector <- contact_data$demography$population

# set demography vector names
names(demography_vector) <- contact_data$demography$age.group

# scale by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)
n_risk_grps <- 1L

# prepare p_susceptibility and susceptibility
p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)

# set susceptibility of group 1 to ZERO
susceptibility <- matrix(
  data = c(0.0, 1.0, 1.0), nrow = n_demo_grps, ncol = n_risk_grps
)

#### Test output using iterative solver ####
# prepare final_size output
epi_outcome <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "iterative"
)

# Test that final_size values are within range and have correct length
test_that("Solver speed-ups return correct final size values", {
  # check for bad numeric, NAN, or infinite values
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

  # expect that first group has p_infected = 0.0
  expect_identical(
    epi_outcome$p_infected[1], 0.0
  )
})

#### Test output using Newton solver ####
# prepare final_size output
epi_outcome <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "newton"
)

# Test that final_size values are within range and have correct length
test_that("Newton solver speed-up returns correct value", {
  # expect that first group has p_infected = 0.0
  expect_identical(
    epi_outcome$p_infected[1], 0.0
  )
})
