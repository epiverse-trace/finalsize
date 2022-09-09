test_that("Check final_size_cpp calculation works", {
  c_matrix <- matrix(0.2, 2, 2)
  d_vector <- c(0.5, 0.5)
  p_suscep <- 1.0
  p_initial_infections <- 0.0015

  epi_final_size <- final_size_cpp(
    r0 = 1.3,
    contact_matrix = c_matrix,
    demography_vector = d_vector,
    prop_initial_infected = p_initial_infections,
    prop_suscep = p_suscep
  )

  # basic tests for output
  testthat::expect_identical(
    length(d_vector), length(epi_final_size)
  )
  # test that output is a numeric vector
  testthat::expect_vector(
    epi_final_size,
    ptype = numeric()
  )
  # test that the final size is plausible
  testthat::expect_true(
    all(epi_final_size <= 1.0)
  )
  # expect that final size is greater than p_initial_infections
  testthat::expect_true(
    all(epi_final_size >= p_initial_infections)
  )
})

test_that("Check final_size_cpp with Polymod data", {
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
  # basic tests for output
  testthat::expect_identical(
    length(d_vector), length(epi_final_size)
  )
  # test that output is a numeric vector
  testthat::expect_vector(
    epi_final_size,
    ptype = numeric()
  )
  # test that the final size is plausible
  testthat::expect_true(
    all(epi_final_size <= 1.0)
  )
  # expect that final size is greater than p_initial_infections
  testthat::expect_true(
    all(epi_final_size >= p_initial_infections)
  )
})
