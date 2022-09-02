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
  
  epi_final_size <- final_size(
    r0 = 2,
    contact_matrix = c_matrix,
    demography_vector = d_vector,
    prop_suscep = p_suscep
  )
  
  # Run final size model
  testthat::expect_identical(
    length(p_suscep), length(epi_final_size)
  )
})
