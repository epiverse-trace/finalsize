# Calculates the upper limit of final size given the r0
# The upper limit is given by a well mixed population
upper_limit <- function(r0) {
  f <- function(par) {
    abs(1 - exp(-r0 * par[1]) - par[1])
  }
  opt <- optim(
    par = c(0.5), fn = f,
    lower = c(0), upper = c(1),
    method = "Brent"
  )
  opt
}

test_that("Final size calculations are correct", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  p_initial_infections <- 0.002
  p_suscep <- c(1, 1, 1)

  # r0 greater than 1
  vector_r0 <- seq(1.01, 3.0, 0.5)

  # check that final sizes > initial sizes
  invisible(
    lapply(
      vector_r0, function(r0) {
        fs <- final_size(
          r0 = r0,
          contact_matrix = c_matrix,
          demography_vector = d_vector,
          prop_initial_infected = p_initial_infections,
          prop_suscep = p_suscep
        )

        testthat::expect_true(
          all(fs >= p_initial_infections)
        )

        testthat::expect_true(
          all(fs <= 1.0)
        )
      }
    )
  )
})
