# basic test that solver returns numerics within range
test_that("Epi spread function works", {
  # prepare some data for the solver
  r0 <- 1.3 # not necessary
  contact_matrix <- matrix(r0 / 200.0, 2, 2)
  demography_vector <- rep(100.0, 2)
  demography_vector <- demography_vector / sum(demography_vector)

  # p_susceptibility and susceptibility
  n_risk_grps <- 3
  psusc <- matrix(1, nrow = 2, ncol = n_risk_grps)
  psusc <- t(apply(psusc, 1, \(x) x / sum(x))) # rows sum to 1.0
  susc <- psusc

  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc
  )

  # check names in epi_spread_data
  expected_names <- c(
    "contact_matrix", "demography_vector", "p_susceptibility",
    "susceptibility"
  )
  expect_equal(
    setdiff(
      names(epi_spread_data), expected_names
    ), character()
  )

  # check that solver returns correct types
  expect_true(
    is.matrix(epi_spread_data[["contact_matrix"]])
  )
  expect_vector(
    epi_spread_data[["demography_vector"]],
    ptype = numeric()
  )
  expect_true(
    is.matrix(epi_spread_data[["p_susceptibility"]])
  )
  expect_equal(
    ncol(epi_spread_data[["p_susceptibility"]]),
    1
  )
  expect_equal(
    nrow(epi_spread_data[["p_susceptibility"]]),
    prod(dim(psusc))
  )
  expect_equal(
    ncol(epi_spread_data[["susceptibility"]]),
    1
  )
  expect_equal(
    nrow(epi_spread_data[["susceptibility"]]),
    prod(dim(susc))
  )
  # expect that replicated demography is always lower than input
  # when there are multiple risk groups
  expect_true(
    all(epi_spread_data[["demography_vector"]] < demography_vector)
  )
})
