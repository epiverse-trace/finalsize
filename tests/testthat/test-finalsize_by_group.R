test_that("Check basic final size function works", {
  # checking epi spread function from finalsize
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40),
    symmetric = TRUE
  )

  demography_vector <- contact_data$demography$proportion

  # Scale contact matrix to demography
  contact_matrix <- apply(
    contact_data$matrix, 1, function(r) r / demography_vector
  )
  # This is needed to pass isSymmetric. Maybe isSymmetric does not work
  # properly with named matrices?
  # TODO: fix the function so that this is not needed any more
  contact_matrix <- matrix(contact_matrix, ncol = 3)
  testthat::expect_true(isSymmetric(contact_matrix))

  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 1, 3)

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  testthat::expect_type(
    epi_outcome, "double"
  )
})

test_that("Check final size calculation is correct in simple case", {
  r0 <- 2
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2) |> as.matrix()
  psusc <- rep(1.0, 2) |> as.matrix()
  susc <- rep(1.0, 2) |> as.matrix()

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  epi_outcome_known <- 1 - exp(-r0 * epi_outcome)

  testthat::expect_equal(
    epi_outcome, epi_outcome_known
  )
})

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
  r0 <- 2.0

  # a p_susceptibility matrix
  p_susc <- matrix(0, nrow(contact_matrix), 4) # four susceptibiliy groups
  # fill p_susceptibility columns manually
  p_susc[, 1] <- 0.7
  p_susc[, 2] <- 0.1
  p_susc[, 3] <- 0.1
  p_susc[, 4] <- 0.1

  # a susceptibility matrix
  susc <- matrix(0, nrow(contact_matrix), 4)
  susc[, 1] <- 1.0
  susc[, 2] <- 0.7
  susc[, 3] <- 0.4
  susc[, 4] <- 0.1

  epi_outcome <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susc,
    susceptibility = susc
  )

  testthat::expect_type(
    epi_outcome, "double"
  )

  testthat::expect_equal(
    length(epi_outcome), length(demography_vector)
  )

  # TO DO: ADD CHECK FOR CORRECT ANSWER
})

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
  testthat::expect_error(
    final_size_grps_cpp(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: contact matrix must have as many rows as demography groups"
  )

  demography_vector <- demography_vector[-length(demography_vector)]
  p_susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and p_susceptibility
  testthat::expect_error(
    final_size_grps_cpp(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: p_susceptibility must have as many rows as demography"
  )

  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 1, nrow = 4)

  # expect error on demography vector and susceptibility
  testthat::expect_error(
    final_size_grps_cpp(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: susceptibility must have as many rows as demography groups"
  )

  p_susceptibility <- matrix(1, ncol = 1, nrow = 3)
  susceptibility <- matrix(1, ncol = 2, nrow = 3)

  # expect error on p_susceptibility and susceptibility
  testthat::expect_error(
    final_size_grps_cpp(
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    regexp = "Error: p_susceptibility and susceptibility must be matrices of the same dims"
  )
})
