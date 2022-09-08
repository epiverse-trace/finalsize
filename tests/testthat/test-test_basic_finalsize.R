test_that("Check basic final size calculation works", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  p_suscep <- c(1, 1, 1)
  p_initial_infections <- 0.0015

  epi_final_size <- final_size(
    r0 = 2,
    contact_matrix = c_matrix,
    demography_vector = d_vector,
    prop_initial_infected = p_initial_infections,
    prop_suscep = p_suscep
  )

  # Run final size model
  testthat::expect_identical(
    length(p_suscep), length(epi_final_size)
  )
})

test_that("Check final size r0 inputs", {
  c_matrix <- matrix(0.5, 2, 2)
  d_vector <- rep(0.5, 2)
  p_suscep <- c(1)
  p_initial_infections <- 0.0015

  # expect error for negative r0
  testthat::expect_error(
    final_size(
      r0 = -1.0,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(Element 1 is not >= 1)"
  )

  # expect error for r0 < 1
  testthat::expect_error(
    final_size(
      r0 = 0.99,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(Element 1 is not >= 1)"
  )

  # expect error for infinite r0
  testthat::expect_error(
    final_size(
      r0 = Inf,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(Must be finite)"
  )
})

test_that("Check failure for unequal contact matrix and demography", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- c(0.5, 0.5) # manually created
  p_suscep <- c(1, 1, 1)
  p_initial_infections <- c(0.1, 0.5, 0.5)

  testthat::expect_error(
    final_size(
      r0 = 2,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(demography vector)*(same size)*(contact matrix)"
  )
})

test_that("Check failure for unequal demography and susceptibility vectors", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- c(0.5, 0.5, 0.5) # manually created
  p_suscep <- c(1, 1, 1, 1) # one extra susceptibility value
  p_initial_infections <- c(0.0015)

  testthat::expect_error(
    final_size(
      r0 = 2,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(demography vector)*(same size)*(susceptibility vector)"
  )
})

test_that("Check message for multiple initial infection proportions", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion # manually created
  p_suscep <- c(1, 1, 1)
  p_initial_infections <- c(0.1, 0.5, 0.5)

  testthat::expect_message(
    final_size(
      r0 = 2,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(different)*(prop_initial_infected)*(age group)"
  )
})

test_that("Check error for unequal proportions initial infections", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion # manually created
  p_suscep <- c(1, 1, 1)
  p_initial_infections <- c(0.1, 0.5, 0.5, 0.1)

  testthat::expect_error(
    final_size(
      r0 = 2,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "(vector)*(prop_initial_infected)*(same size)*(demography vector)"
  )
})
