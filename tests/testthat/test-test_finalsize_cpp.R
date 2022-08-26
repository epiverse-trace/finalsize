test_that("Check basic final size calculation works", {
  # checking epi spread function from finalsize
  polymod <- socialmixr::polymod
  contact_matrix <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40)
  )
  demography <- contact_matrix$participants$proportion
  p_susceptibility <- matrix(1, 3, 3)
  susceptibility = matrix(1, 3, 3)

  epi_outcome = final_size_cpp(
    contact_matrix = contact_matrix,
    demography = demography,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  testthat::expect_type(
    epi_outcome, "double"
  )

})
