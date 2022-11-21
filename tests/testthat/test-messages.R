#### Test that errors, warnings, and messages are triggered correctly ####

# Prepare common elements for testing
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

# scale by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- contact_matrix / demography_vector

p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
susceptibility <- matrix(1, ncol = 1, 3)


#### Check for demography vector length mismatch ####
test_that("Check for contact_matrix -- demography vector mismatch error", {

  # wrong demography vector
  demography_vector <- c(demography_vector, 100)

  # expect error on demography vector and contact matrix
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: contact matrix must have as many rows as demography groups"
  )
})

test_that("Check for p_susceptibility -- demography vector mismatch error", {
  p_susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and p_susceptibility
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp = "Error: p_susceptibility must have as many rows as demography"
  )
})

test_that("Check for susceptibility -- demography mismatch error", {
  susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and susceptibility
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp = "Error: susceptibility must have as many rows as demography groups"
  )
})

test_that("Check for susceptibility -- p_susceptibility mismatch error", {
  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 2, nrow = 3)

  # expect error on p_susceptibility and susceptibility
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp =
      "Error: susceptibility must have same dimensions as p_susceptibility"
  )
})

test_that("Check for solver option error", {
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "some wrong solver option"
    ),
    regexp =
      "(Error)*(should be one of)*(iterative)*(newton)"
  )
})

test_that("Check for warning when error is much larger than tolerance", {
  # for the iterative solver
  expect_warning(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list(
        iterations = 2,
        tolerance = 1e-12
      )
    ),
    regexp = "Solver error > 100x solver tolerance, try increasing iterations"
  )

  # for the newton solver
  expect_warning(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "newton",
      control = list(
        iterations = 2,
        tolerance = 1e-12
      )
    ),
    regexp =
      "Solver error > 100x solver tolerance, try increasing iterations"
  )
})

test_that("Check that arguments are of the correct type", {
  # expect errors when wrong argument types are passed
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = as.vector(contact_matrix),
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: contact matrix must be a matrix"
  )
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = as.matrix(demography_vector),
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: demography vector must be a numeric vector"
  )
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = as.vector(p_susceptibility),
      susceptibility = susceptibility
    ),
    regexp = "Error: p_susceptibility must be a matrix"
  )
  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = as.vector(susceptibility),
      p_susceptibility = p_susceptibility
    ),
    regexp = "Error: susceptibility must be a matrix"
  )
})

# Check the contents of the control list
test_that("Check that eigenvalue checking works", {
  r0 <- 2
  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 1, 3)

  expect_error(
    final_size(
      r0 = r0,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      control = list(some_other_name = 10000)
    ),
    regexp = "Error: control list names can only be:"
  )
})
