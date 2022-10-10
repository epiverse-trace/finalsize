#' @title Prepare contact and demography data for multiple risk groups
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param p_susceptibility A matrix giving the probability that an individual
#' in demography group $i$ is in risk (or susceptibility) group $j$.
#' Each row represents the overall distribution of individuals in demographic
#' group $i$ across risk groups, and each row *must sum to 1.0*.
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group $i$ and risk group $j$.
#'
#' @return A list object with named elements: "contact_matrix",
#' "demography_vector", "p_susceptibility_", and "susceptibility".
#' The contact matrix is replicated row and column wise for each risk group
#' and the demography vector is replicated for each risk group.
epi_spread <- function(contact_matrix,
                       demography_vector,
                       p_susceptibility,
                       susceptibility) {
  # count risk groups
  n_susc_groups <- ncol(p_susceptibility)
  # make lps, a 1 col matrix of all p_susc values
  lps <- as.vector(p_susceptibility)
  # replicate the demography vector and multiply by p_susceptibility
  demography_vector_spread <- rep(demography_vector, n_susc_groups)
  demography_vector_spread <- demography_vector_spread * lps

  # replicate contact matrix
  contact_matrix_spread <- kronecker(
    X = matrix(1, nrow = n_susc_groups, ncol = n_susc_groups),
    Y = contact_matrix
  )

  # unroll the susceptibility matrix
  susceptibility_ <- as.vector(susceptibility)

  list(
    contact_matrix = contact_matrix_spread,
    demography_vector = demography_vector_spread,
    susceptibility = susceptibility_
  )
}
