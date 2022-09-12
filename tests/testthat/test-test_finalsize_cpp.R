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

  epi_final_size <- final_size_cpp(
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

test_that("Check final_size_cpp r0 values increase final size", {
  c_matrix <- matrix(0.2, 2, 2)
  d_vector <- c(0.5, 0.5)
  p_suscep <- c(1.0, 1.0)
  p_initial_infections <- c(0.0015)

  r0_high <- 3.0
  r0_low <- 1.1

  # basic tests for output
  final_sizes <- lapply(
    list(r0_high, r0_low),
    function(r0_) {
      final_size_cpp(
        r0 = r0_,
        contact_matrix = c_matrix,
        demography_vector = d_vector,
        prop_initial_infected = p_initial_infections,
        prop_suscep = p_suscep
      )
    }
  )

  # test that final sizes for higher r0 are higher than low r0
  testthat::expect_true(
    all(final_sizes[[1]] > final_sizes[[2]])
  )
})

test_that("Check final_size_cpp calculation works", {
  c_matrix <- matrix(0.2, 2, 2)
  d_vector <- c(0.5, 0.5)
  p_suscep <- c(1.0, 1.0, 1.0)
  p_initial_infections <- c(0.0015)

  p_initial_infections_wrong <- c(0.01, 0.01, 0.01)
  p_suscep_correct <- 0.01

  # basic tests for output
  testthat::expect_error(
    final_size_cpp(
      r0 = 1.3,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections,
      prop_suscep = p_suscep
    ),
    regexp = "Error: prop_suscep must be same size as demography vector"
  )

  # basic tests for output
  testthat::expect_error(
    final_size_cpp(
      r0 = 1.3,
      contact_matrix = c_matrix,
      demography_vector = d_vector,
      prop_initial_infected = p_initial_infections_wrong,
      prop_suscep = p_suscep_correct
    ),
    regexp = "Error: prop_initial_infection must be same size as demography"
  )
})
