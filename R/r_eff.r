#' @title Effective \eqn{R_0} in a heterogeneous population
#'
#' @description `r_eff` calculates the effective reproductive number \eqn{R_eff}
#' in a population with heterogeneous mixing, and with heterogeneous
#' susceptibility to infection such as due to immunisation.
#'
#'
#' @param r0 The basic reproductive number \eqn{R_0} of the infection.
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i} (model will normalise
#' if needed).
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group \eqn{i} and risk group \eqn{j}.
#' @param p_susceptibility A matrix giving the probability that an individual
#' in demography group \eqn{i} is in risk (or susceptibility) group \eqn{j}.
#' Each row represents the overall distribution of individuals in demographic
#' group \eqn{i} across risk groups, and each row *must sum to 1.0*.
#'
#' @keywords R0
#' @return A single number of the effective reproductive number of the infection
#' in the population.
#' @export
#' @examples
#' # load example POLYMOD data included in the package
#' data(polymod_uk)
#' r0 <- 2.0
#' contact_matrix <- polymod_uk$contact_matrix
#' demography_vector <- polymod_uk$demography_vector
#'
#' # define the number of age and susceptibility groups
#' n_demo_grps <- length(demography_vector)
#' n_risk_grps <- 3
#'
#' # In this example, all risk groups from all age groups are fully
#' # susceptible
#' susceptibility <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#'
#' p_susceptibility <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#' # p_susceptibility rows must sum to 1.0
#' p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
#'
#' # calculate R_effective
#' r_eff(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility
#' )
r_eff <- function(r0,
                  contact_matrix,
                  demography_vector,
                  susceptibility,
                  p_susceptibility) {
  # check arguments input
  stopifnot(
    "Error: contact matrix must be a matrix" =
      (is.matrix(contact_matrix)),
    "Error: demography vector must be a numeric vector" =
      (is.vector(demography_vector) & is.numeric(demography_vector)),
    "Error: p_susceptibility must be a matrix" =
      (is.matrix(p_susceptibility)),
    "Error: susceptibility must be a matrix" =
      (is.matrix(susceptibility)),
    "Error: contact matrix must have as many rows as demography groups" =
      (nrow(contact_matrix) == length(demography_vector)),
    "Error: p_susceptibility must have as many rows as demography groups" =
      (nrow(p_susceptibility) == length(demography_vector)),
    "Error: susceptibility must have as many rows as demography groups" =
      (nrow(susceptibility) == length(demography_vector)),
    "Error: susceptibility must have same dimensions as p_susceptibility" =
      (all(dim(p_susceptibility) == dim(susceptibility))),
    "Error: p_susceptibility rows must sum to 1.0" =
      (
        all(abs(rowSums(p_susceptibility) - 1) < 1e-6)
      ),
    "Error: contact matrix must have a maximum real eigenvalue of 1.0" =
      (
        abs(max(Re(eigen(contact_matrix * demography_vector)$values) - 1.0)) <
          1e-6
      )
  )

  # count risk groups
  n_susc_groups <- ncol(p_susceptibility)
  # make lps, a vector of all p_susc values
  lps <- as.vector(p_susceptibility)

  # replicate the demography vector and multiply by p_susceptibility
  demography_vector_spread <- rep(demography_vector, n_susc_groups)
  demography_vector_spread <- demography_vector_spread * lps

  # replicate contact matrix
  contact_matrix_spread <- kronecker(
    X = matrix(1, nrow = n_susc_groups, ncol = n_susc_groups),
    Y = contact_matrix * r0
  )

  # set contact_matrix values to zero if there are no contacts among ...
  # demography groups, or if demography groups are empty
  i_here <- demography_vector_spread == 0 | as.vector(susceptibility) == 0 |
    rowSums(contact_matrix_spread) == 0
  contact_matrix_spread[i_here, i_here] <- 0.0

  # process contact matrix in solver specific ways
  contact_matrix_spread <- t(
    t(contact_matrix_spread * as.vector(susceptibility)) *
      demography_vector_spread
  )

  # get effective R
  r_eff_ <- max(Re(eigen(contact_matrix_spread)$values))

  # return r_eff_
  r_eff_
}
