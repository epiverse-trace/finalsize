# Check final_size returns data frame with correct group names

# Prepare data
r0 <- 2.0
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  split = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# set demography vector names
names(demography_vector) <- contact_data$demography$age.group

# Scale contact matrix by eigenvalue and demography to demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- apply(
  contact_matrix, 1, function(r) r / demography_vector
)

n_demo_grps <- length(demography_vector)
n_risk_grps <- 2

# prepare p_susceptibility and susceptibility
psusc <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
psusc <- psusc / rowSums(psusc)
susc <- matrix(
  data = c(1.0, 0.5), nrow = n_demo_grps, ncol = n_risk_grps,
  byrow = TRUE
)
colnames(susc) <- c("susceptible", "immunised")

# Test for demography names from named demography vector
test_that("Finalsize returns correct demography names", {
  epi_outcome <- final_size(
    contact_matrix = r0 * contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative"
  )

  expect_s3_class(
    epi_outcome, "data.frame"
  )

  # check for names
  expect_equal(
    unique(epi_outcome$demo_grp),
    names(demography_vector)
  )
  expect_equal(
    unique(epi_outcome$susc_grp), c("susceptible", "immunised")
  )
})

# Test for demographic names from contact matrix
test_that("Finalsize takes demography names from contact data", {
  rownames(contact_matrix) <- names(demography_vector)
  names(demography_vector) <- NULL
  colnames(susc) <- c("susceptible", "immunised")

  epi_outcome <- final_size(
    contact_matrix = r0 * contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative"
  )

  expect_s3_class(
    epi_outcome, "data.frame"
  )

  # check for names
  expect_equal(
    unique(epi_outcome$demo_grp),
    rownames(contact_matrix)
  )
  expect_equal(
    unique(epi_outcome$susc_grp), c("susceptible", "immunised")
  )
})

# Test for demographic names from contact matrix
test_that("Finalsize generates group names", {
  names(demography_vector) <- NULL
  rownames(contact_matrix) <- NULL
  colnames(susc) <- NULL

  epi_outcome <- final_size(
    contact_matrix = r0 * contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = psusc,
    susceptibility = susc,
    solver = "iterative"
  )

  expect_s3_class(
    epi_outcome, "data.frame"
  )

  # check for names
  expect_equal(
    unique(epi_outcome$demo_grp),
    sprintf("demo_grp_%i", seq(length(demography_vector)))
  )
  expect_equal(
    unique(epi_outcome$susc_grp),
    sprintf(
      "susc_grp_%i", seq_len(ncol(susc))
    )
  )
})
