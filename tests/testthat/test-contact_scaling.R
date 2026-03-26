test_that("Contact scaling works in `final_size()`", {
  r0 <- 1.5

  # simple case
  expect_no_condition(
    final_size(r0, contact_scaling = 0.5)
  )

  scaling <- 0.95
  fs <- final_size(r0)
  fs_scaled <- final_size(r0, contact_scaling = scaling)

  expect_lt(
    fs_scaled$p_infected, fs$p_infected
  )
  expect_false(
    identical(fs_scaled$p_infected, fs$p_infected * scaling)
  )

  expect_error(
    final_size(r0, contact_scaling = -1.0)
  )
  expect_error(
    final_size(r0, contact_scaling = 1.1)
  )
  expect_error(
    final_size(r0, contact_scaling = c(1, 1))
  )

  # complex case
  demography <- rep(1000, 4)
  cm <- (matrix(1, 4, 4) / 4) / demography
  p_susc <- matrix(1, 4, 1)
  susc <- matrix(1, 4, 1)
  scaling <- rep(0.95, 4)

  expect_no_condition(
    final_size(r0, cm, demography, susc, p_susc, scaling)
  )

  expect_true(
    all(final_size(r0, cm, demography, susc, p_susc, scaling)$p_infected <
      final_size(r0, cm, demography, susc, p_susc)$p_infected)
  )
})

test_that("Contact scaling works in `r_eff()`", {
  r0 <- 1.5

  # simple case
  expect_no_condition(
    r_eff(r0, matrix(1), 1, matrix(1), matrix(1), contact_scaling = 0.5)
  )

  scaling <- 0.95
  reff_scaled <- r_eff(r0, matrix(1), 1, matrix(1), matrix(1), scaling)

  expect_lt(
    reff_scaled, r0
  )

  expect_error(
    r_eff(r0, matrix(1), 1, matrix(1), matrix(1), contact_scaling = -1.0)
  )
  expect_error(
    r_eff(r0, matrix(1), 1, matrix(1), matrix(1), contact_scaling = 1.1)
  )
  expect_error(
    r_eff(r0, matrix(1), 1, matrix(1), matrix(1), contact_scaling = c(1, 1))
  )

  # complex case
  demography <- rep(1000, 4)
  cm <- (matrix(1, 4, 4) / 4) / demography
  p_susc <- matrix(1, 4, 1)
  susc <- matrix(1, 4, 1)
  scaling <- rep(0.95, 4)

  expect_no_condition(
    final_size(r0, cm, demography, susc, p_susc, scaling)
  )

  expect_true(
    all(final_size(r0, cm, demography, susc, p_susc, scaling)$p_infected <
      final_size(r0, cm, demography, susc, p_susc)$p_infected)
  )
})
