#### Tests for the conversion from R0 to lambda and vice versa ####
# Get example dataset and prepare contact matrix and demography
data(polymod_uk)
contact_matrix <- polymod_uk$contact_matrix
demography_vector <- polymod_uk$demography_vector

# define the number of age and susceptibility groups
n_demo_grps <- length(demography_vector)
n_risk_grps <- 3

# In this example, all risk groups from all age groups are fully
# susceptible
susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)

p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
# p_susceptibility rows must sum to 1.0
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)

# define infectious period of 5 days
infectious_period <- 5

# tests for lambda to R0
test_that("Lambda is converted to R0", {
  # basic test for conversion
  lambda <- 0.3

  r0 <- lambda_to_r0(
    lambda = lambda, contact_matrix = contact_matrix,
    susceptibility = susceptibility, p_susceptibility = p_susceptibility,
    demography_vector = demography_vector, infectious_period = infectious_period
  )
  expect_type(r0, "double")
  expect_length(r0, 1)

  # tests for correctness - increasing lambda should increase r0
  r0_vec <- vapply(
    c(0.3, 0.5), lambda_to_r0,
    FUN.VALUE = numeric(1),
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    infectious_period = infectious_period
  )
  expect_gt(
    r0_vec[2], r0_vec[1]
  )

  # expect that heterogeneous susceptibility leads to lower r0
  susc_variable <- matrix(
    data = c(0.1, 0.5, 1.0),
    nrow = n_demo_grps, ncol = n_risk_grps,
    byrow = TRUE
  )
  r0_variable_susc <- lambda_to_r0(
    lambda = lambda, contact_matrix = contact_matrix,
    susceptibility = susc_variable,
    p_susceptibility = p_susceptibility,
    demography_vector = demography_vector, infectious_period = infectious_period
  )
  expect_lt(
    r0_variable_susc, r0
  )
})

# tests for R0 to lambda
test_that("R0 is converted to lambda", {
  # basic test for conversion
  r0 <- 1.5

  lambda <- r0_to_lambda(
    r0 = r0,
    contact_matrix = contact_matrix, demography_vector = demography_vector,
    susceptibility = susceptibility, p_susceptibility = p_susceptibility,
    infectious_period = infectious_period
  )
  expect_type(lambda, "double")
  expect_length(lambda, 1)

  # tests for correctness - increasing lambda should increase r0
  lambda_vec <- vapply(
    c(1.3, 1.5), r0_to_lambda,
    FUN.VALUE = numeric(1),
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    infectious_period = infectious_period
  )
  expect_gt(
    lambda_vec[2], lambda_vec[1]
  )

  # expect that heterogeneous susceptibility leads to lower lambda
  susc_variable <- matrix(
    data = c(0.9, 1.0, 1.0),
    nrow = n_demo_grps, ncol = n_risk_grps,
    byrow = TRUE
  )
  lambda_variable_susc <-
    r0_to_lambda(
      r0 = r0, contact_matrix = contact_matrix,
      susceptibility = susc_variable,
      p_susceptibility = p_susceptibility,
      demography_vector = demography_vector,
      infectious_period = infectious_period
    )
  expect_lt(
    lambda_variable_susc, lambda
  )
})

# tests for equivalence
test_that("Conversion functions are inverses of each other", {
  # convert lambda to R0 and back
  lambda <- 0.3
  r0 <- lambda_to_r0(
    lambda, contact_matrix, demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    infectious_period
  )
  expect_identical(
    r0_to_lambda(
      r0 = r0,
      contact_matrix, demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      infectious_period
    ),
    lambda
  )
})
