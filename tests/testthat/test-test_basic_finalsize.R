# Calculates the upper limit of final size given the r0
# The upper limit is given by a well mixed population
upper_limit <- function(r0) {
  f <- function(par) {
    abs(1 - exp(-r0 * par[1]) - par[1])
  }
  optim(par = c(0.5), fn = f, lower = c(0), upper = c(1), 
    method = "Brent") -> opt
  opt
}

test_that("Check basic final size calculation works", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  suscep <- c(1, 1, 1)
  r0_value <- 2.0

  epi_final_size <- final_size(
    r0 = r0_value,
    contact_matrix = c_matrix,
    demography = d_vector,
    susceptibility = suscep
  )

  # Run final size model
  testthat::expect_identical(
    length(suscep), length(epi_final_size)
  )

  # This should be valid for all values of r0
  testthat::expect_true(all(epi_final_size < 1))
  testthat::expect_true(all(epi_final_size > 0))
  # Calculate the final size in the population
  pi <- sum(epi_final_size*d_vector)/sum(d_vector)
  max_pi <- upper_limit(r0_value)
  testthat::expect_equal(max_pi$convergence, 0)
  testthat::expect_lte(pi, max_pi$par)

  # Lower realistic limit given an r0 of 2
  testthat::expect_true(all(epi_final_size > 0.7))

  # Check that lower r0 values result in lower final size values
  r0_value_low <- 1.5

  epi_final_size_low <- final_size(
    r0 = r0_value_low,
    contact_matrix = c_matrix,
    demography = d_vector,
    susceptibility = p_suscep
  )
  testthat::expect_true(all(epi_final_size_low  < epi_final_size))
  testthat::expect_true(all(epi_final_size_low >= 0))
})

test_that("Check that R0 < 1 results in a final size of zero", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  suscep <- c(1, 1, 1)
  r0_value <- 0.9

  epi_final_size <- final_size(
    r0 = r0_value,
    contact_matrix = c_matrix,
    demography = d_vector,
    susceptibility = suscep
  )

  testthat::expect_identical(
    length(suscep), length(epi_final_size)
  )

  testthat::expect_true(all(epi_final_size == 0))
})

test_that("Check that final size works properly if one age group is not susceptible", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )

  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  suscep <- c(1, 1, 0)
  r0_value <- 2.9

  epi_final_size <- final_size(
    r0 = r0_value,
    contact_matrix = c_matrix,
    demography = d_vector,
    susceptibility = suscep
  )

  testthat::expect_identical(
    length(suscep), length(epi_final_size)
  )

  testthat::expect_true(all(epi_final_size[1:2] > 0))
  testthat::expect_true(epi_final_size[3] == 0)

  pi <- sum(epi_final_size*d_vector)/sum(d_vector)
  max_pi <- upper_limit(r0_value)
  testthat::expect_equal(max_pi$convergence, 0)
  testthat::expect_lte(pi, max_pi$par)
})
