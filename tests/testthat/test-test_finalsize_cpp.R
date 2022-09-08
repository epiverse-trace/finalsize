test_that("Check final_size_cpp calculation works", {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  c_matrix <- t(contact_data$matrix)
  d_vector <- contact_data$participants$proportion
  p_suscep <- c(1)
  p_initial_infections <- 0.0015

  for(i in seq(ncol(c_matrix))) {
    c_matrix[, i] = c_matrix[, i] / d_vector[i]
  }

  epi_final_size <- final_size_cpp(
    r0 = 1.3,
    contact_matrix = c_matrix,
    demography_vector = d_vector,
    prop_initial_infected = p_initial_infections,
    prop_suscep = p_suscep
  )

  # Run final size model
  testthat::expect_identical(
    length(d_vector), length(epi_final_size)
  )
})
