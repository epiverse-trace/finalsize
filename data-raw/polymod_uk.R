## code to prepare `polymod_uk` dataset goes here

# load POLYMOD data for the UK from {socialmixr}
polymod_uk <- socialmixr::polymod
polymod_uk <- socialmixr::contact_matrix(
  polymod_uk,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)

# get contact matrix and demography vector
contact_matrix <- t(polymod_uk$matrix)
demography_vector <- polymod_uk$demography$population

# scale contacts by maximum real eigenvalue and divide by demography
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)
contact_matrix <- contact_matrix / demography_vector

polymod_uk <- list(
  contact_matrix = contact_matrix,
  demography_vector = demography_vector
)

usethis::use_data(polymod_uk, overwrite = TRUE)
