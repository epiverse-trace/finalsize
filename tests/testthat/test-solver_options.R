#### Tests for solver options ####

#### Prepare data ####
r0 <- 2.0
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# set demography vector names
names(demography_vector) <- contact_data$demography$age.group

# scale by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)
n_risk_grps <- 2

# prepare p_susceptibility and susceptibility
p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
susceptibility <- matrix(
  data = c(1.0, 0.5), nrow = n_demo_grps, ncol = n_risk_grps,
  byrow = TRUE
)
colnames(susceptibility) <- c("susceptible", "immunised")

#### Test output using iterative solver options ####
# prepare final_size output
epi_outcome <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "iterative",
  control = list(
    adapt_step = FALSE # default option is TRUE
  )
)

# Test that final_size values are within range and have correct length
test_that("Iterative solver with adaptive steps turned off works", {
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
})
