# Check final_size works with Newton solver
# check for errors and messages
test_that("Check for errors and messages", {
  # checking epi spread function from finalsize
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40),
    symmetric = TRUE
  )
  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 1, 3)

  demography_vector <- contact_data$demography$proportion

  # 'wrong' demography vector
  demography_vector <- c(demography_vector, 100)

  contact_matrix <- matrix(contact_data$matrix, ncol = 3)

  # expect error on demography vector and contact matrix
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp = "Error: contact matrix must have as many rows as demography groups"
  )

  demography_vector <- demography_vector[-length(demography_vector)]
  p_susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and p_susceptibility
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp = "Error: p_susceptibility must have as many rows as demography"
  )

  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and susceptibility
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp = "Error: susceptibility must have as many rows as demography groups"
  )

  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 2, nrow = 3)

  # expect error on p_susceptibility and susceptibility
  expect_error(
    final_size(
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

  # expect error when p_susceptibility sum > 1
  p_susceptibility <- matrix(1, ncol = 2, nrow = 3)
  susceptibility <- matrix(1, ncol = 2, nrow = 3)
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list()
    ),
    regexp =
      "Error: p_susceptibility rows must sum to 1.0"
  )

  # expect error when incorrect solver option is passed
  p_susceptibility <- matrix(1, ncol = 2, nrow = 3)
  p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "some wrong solver option",
      control = list()
    ),
    regexp =
      "(Error)*(should be one of)*(iterative)*(newton)"
  )

  # check for warning when error is much larger than tolerance, iterative
  expect_warning(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility,
      solver = "iterative",
      control = list(
        iterations = 2,
        tolerance = 1e-12,
        adapt_step = TRUE,
        step_rate = 1.9
      )
    ),
    regexp =
      "Solver error > 100x solver tolerance, try increasing iterations"
  )

  # check for warning when error is much larger than tolerance, newton
  expect_warning(
    final_size(
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

  # expect errors when wrong argument types are passed
  expect_error(
    final_size(
      contact_matrix = as.vector(contact_matrix),
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: contact matrix must be a matrix"
  )
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = as.matrix(demography_vector),
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: demography vector must be a numeric vector"
  )
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = as.vector(p_susceptibility),
      susceptibility = susceptibility
    ),
    regexp = "Error: p_susceptibility must be a matrix"
  )
  expect_error(
    final_size(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = as.vector(susceptibility),
      p_susceptibility = p_susceptibility
    ),
    regexp = "Error: susceptibility must be a matrix"
  )
})
