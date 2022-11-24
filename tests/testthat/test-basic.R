#### Tests for basic finalsize functionality ####

#### Prepare data ####
r0 <- 2.0
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# set demography vector names
names(demography_vector) <- contact_data$demography$age.group

# scale by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)
n_risk_grps <- 2

# prepare p_susceptibility and susceptibility
p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
susceptibility <- matrix(
  data = c(1.0, 0.5), nrow = n_demo_grps, ncol = n_risk_grps,
  byrow = TRUE
)
colnames(susceptibility) <- c("susceptible", "immunised")

# prepare final_size output
epi_outcome <- final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  p_susceptibility = p_susceptibility,
  susceptibility = susceptibility,
  solver = "iterative"
)

#### Test that final_size returns a data.frame with correct columns ####
test_that("Finalsize returns correct demography names", {
  # expect data.frame as output
  expect_s3_class(
    epi_outcome, "data.frame"
  )

  # expect that data.frame has correct names
  expect_identical(
    colnames(epi_outcome),
    c("demo_grp", "susc_grp", "susceptibility", "p_infected")
  )

  # check for names
  expect_named(
    demography_vector,
    unique(epi_outcome$demo_grp)
  )
  expect_setequal(
    epi_outcome$susc_grp, c("susceptible", "immunised")
  )

  # check for number of rows of data.frame
  expect_identical(
    nrow(epi_outcome),
    as.integer(n_demo_grps * n_risk_grps)
  )
})

# Test that final_size values are within range and have correct length
test_that("final_size output is correct length and within range 0 - 1", {
  # check for bad numeric, NAN, or infinite values
  # check that values are not NaN
  expect_false(
    any(is.nan(epi_outcome$p_infected))
  )
  # check that solver returns no nas
  expect_false(
    anyNA(epi_outcome$p_infected)
  )
  # check that solver returns no inf
  expect_false(
    any(is.infinite(epi_outcome$p_infected))
  )
  # check that solver returns values within range
  expect_true(
    all(epi_outcome$p_infected >= 0.0)
  )
  expect_true(
    all(epi_outcome$p_infected <= 1.0)
  )
  # check for number of demography groups returned
  expect_identical(
    n_demo_grps,
    length(unique(epi_outcome$demo_grp))
  )
})

#### Test that final_size returns default group names when required ####

test_that("final_size returns default group names", {
  names(demography_vector) <- NULL
  colnames(contact_matrix) <- NULL
  rownames(contact_matrix) <- NULL
  colnames(susceptibility) <- NULL

  # prepare final_size output with default names
  epi_outcome <- final_size(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility,
    solver = "iterative"
  )
  # demography group names are correct
  expect_identical(
    unique(epi_outcome$demo_grp),
    sprintf("demo_grp_%i", seq_len(n_demo_grps))
  )

  # susceptibility group names are correct
  expect_identical(
    unique(epi_outcome$susc_grp),
    sprintf("susc_grp_%i", seq_len(n_risk_grps))
  )
})
