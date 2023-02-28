#### Tests for calculation of effective R0 ####

# load example POLYMOD data included in the package
data(polymod_uk)
r0 <- 2.0
contact_matrix <- polymod_uk$contact_matrix
demography_vector <- polymod_uk$demography_vector

# define the number of age and susceptibility groups
n_demo_grps <- length(demography_vector)
n_risk_grps <- 3

# Initially all individuals fully susceptible
susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
p_susceptibility <- matrix(
  data = 1, nrow = n_demo_grps, ncol = n_risk_grps
)
# p_susceptibility rows must sum to 1.0
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)

# run basic tests
test_that("`r_eff` returns correct values", {
  # expect that r_eff() returns the same r0 as passed when all are susceptible
  r_eff_all_susc <- r_eff(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility
  )
  expect_identical(r_eff_all_susc, r0, tolerance = 1e-5)

  # expect that r_eff() returns a lower value when susceptibility is reduced
  susceptibility <- matrix(
    data = 0.9, nrow = n_demo_grps, ncol = n_risk_grps
  )
  r_eff_susc_lower <- r_eff(
    r0 = r0,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility
  )
  expect_lt(r_eff_susc_lower, r0)
})
