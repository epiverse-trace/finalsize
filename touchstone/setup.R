library(finalsize)
library("socialmixr")
library("withr")

#### Set up population characteristics ####
# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = seq(0, 80, 10), # final age bin is 70+
  symmetric = TRUE
)

# get the contact matrix and demography data
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# scale the contact matrix so the largest eigenvalue is 1.0
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

# divide each row of the contact matrix by the corresponding demography
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

# get 1000 R0 samples
r0_mean <- 1.5
r0_sd <- 0.01
samples <- 1000
r0_samples <- withr::with_seed(1, rnorm(1000, r0_mean, r0_sd))

# define susceptibility values for age groups, greater for older ages
susc_variable <- matrix(
  data = seq_along(demography_vector) / length(demography_vector)
)
susceptibility <- cbind(
  susc_variable, susc_variable * 0.8, susc_variable * 0.6,
  susc_variable * 0.5
)

# define proportion of each demographic group in each susceptibility group
# assume uniform distribution
n_demo_grps <- nrow(contact_matrix)
p_susceptibility <- matrix(1.0, nrow = n_demo_grps, ncol = ncol(susceptibility))
p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
