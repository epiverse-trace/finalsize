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
  # make p_susceptibility matrix of ones
  p_susceptibility_ <- matrix(
    1.0,
    nrow = prod(dim(p_susceptibility)), ncol = 1
  )
  # make lps, a 1 col matrix of all p_susc values
  lps <- matrix(p_susceptibility, nrow = prod(dim(p_susceptibility)), ncol = 1)
  # replicate the demography vector and multiply by p_susceptibility
  demography_vector_spread <- rep(demography_vector, n_susc_groups)
  demography_vector_spread <- demography_vector_spread * as.vector(lps)

  # replicate contact matrix
  contact_matrix_ <- do.call("rbind", rep(list(contact_matrix), n_susc_groups))
  contact_matrix_ <- do.call(
    "cbind", rep(list(contact_matrix_), n_susc_groups)
  )

  # unroll the susceptibility matrix
  susceptibility_ <- matrix(
    data = as.vector(susceptibility),
    nrow = prod(dim(susceptibility)), ncol = 1
  )

  list(
    contact_matrix = contact_matrix_,
    demography_vector = demography_vector_spread,
    p_susceptibility = p_susceptibility_,
    susceptibility = susceptibility_
  )
}
