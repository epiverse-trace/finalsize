#' Final size of epidemic by susceptibility groups
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param p_susceptibility WIP - what is the correct definition.
#' @param susceptibility WIP - what is the correct definition.
#' @param iterations Number of solver iterations
#' @param tolerance Solver error tolerance.
#' 
#' @return A vector of final sizes by demography group.
#' @export
#' @example 
final_size_grps <- function(contact_matrix,
                                               demography_vector,
                                               p_susceptibility,
                                               susceptibility,
                                               iterations = 1000,
                                               tolerance = 1e-6) {
  # check arguments input
  stopifnot(
    "Error: contact matrix must have as many rows as demography groups" =
      (nrow(contact_matrix) == length(demography_vector))
  )
  stopifnot(
    "Error: p_susceptibility must have as many rows as demography groups" =
      (nrow(p_susceptibility) == length(demography_vector))
  )
  stopifnot(
    "Error: susceptibility must have as many rows as demography groups" =
      (nrow(susceptibility) == length(demography_vector))
  )
  stopifnot(
    "Error: susceptibility must have same dimensions as p_susceptibility" =
      (all(dim(p_susceptibility) == dim(susceptibility)))
  )

  # prepare epi spread object
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  # solve for final size
  pi <- solve_final_size_newton(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]],
    iterations = iterations,
    tolerance = tolerance
  )

  pi_2 <- matrix(pi,
    nrow = length(demography_vector),
    ncol = ncol(p_susceptibility)
  )

  lps <- matrix(p_susceptibility, nrow = prod(dim(p_susceptibility)), ncol = 1)
  v <- as.vector(pi_2) * lps

  pi_3 <- matrix(
    v,
    nrow = length(demography_vector), ncol = ncol(p_susceptibility)
  )

  pi_3
}
