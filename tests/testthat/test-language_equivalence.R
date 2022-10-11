# Check that both R and Cpp implementations give the same answer
test_that("Check language result equivalence", {
  # prepare arguments
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
  r0 <- 1.3
  contact_matrix <- r0 * contact_matrix

  # a p_susceptibility matrix
  p_susceptibility <- matrix(1, nrow(contact_matrix), 1)
  susceptibility <- p_susceptibility

  finalsize_iterative_r <- final_size_grps(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "iterative",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )
  finalsize_newton_r <- final_size_grps(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "newton",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )
  finalsize_iterative_cpp <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "iterative",
    control = list(
      iterations = 10000,
      tolerance = 1e-6,
      step_rate = 1.9,
      adapt_step = TRUE
    )
  )
  finalsize_newton_cpp <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "newton",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )

  # check equivalence between newton solvers R and C++
  expect_equal(
    finalsize_newton_r, finalsize_newton_cpp,
    tolerance = 1e-5
  )
  # check equivalence between iterative solvers R and C++
  expect_equal(
    finalsize_iterative_r, finalsize_iterative_cpp,
    tolerance = 1e-5
  )
})

# Check language equivalence with multiple risk groups
test_that("Check language result equivalence", {
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
  r0 <- 1.3
  contact_matrix <- r0 * contact_matrix

  # a p_susceptibility matrix
  p_susceptibility <- matrix(0.25, nrow(contact_matrix), 4)
  susceptibility <- matrix(c(0.1, 0.7, 0.3, 0.2), nrow(contact_matrix), 4, 
  byrow = T)

  finalsize_iterative_r <- final_size_grps(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "iterative",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )
  finalsize_newton_r <- final_size_grps(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "newton",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )
  finalsize_iterative_cpp <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "iterative",
    control = list(
      iterations = 10000,
      tolerance = 1e-6,
      step_rate = 1.9,
      adapt_step = TRUE
    )
  )
  finalsize_newton_cpp <- final_size_grps_cpp(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility,
    solver = "newton",
    control = list(
      iterations = 10000,
      tolerance = 1e-6
    )
  )

  # check equivalence between newton solvers R and C++
  expect_equal(
    finalsize_newton_r, finalsize_newton_cpp,
    tolerance = 1e-5
  )
  # check equivalence between iterative solvers R and C++
  expect_equal(
    finalsize_iterative_r, finalsize_iterative_cpp,
    tolerance = 1e-5
  )
})
