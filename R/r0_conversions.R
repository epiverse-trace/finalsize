#' @title Calculate \eqn{R_0} from transmission rate (\eqn{\lambda})
#'
#' @description Uses the transmission rate (\eqn{\lambda}), social contacts
#' matrix (\eqn{c}), demography (\eqn{N}), the distribution \eqn{P} of each
#' demographic group \eqn{i} into susceptibility groups \eqn{S}, and the
#' infectious period (\eqn{\gamma}) to calculate the \eqn{R_0} using the
#' following equation.
#' \deqn{R_0 = {Max}(EV(C)) \times \lambda \gamma}
#' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the social contacts matrix scaled by the number of
#' individuals in each demographic and susceptibility group in the population.
#'
#' @param lambda The transmission rate of the disease, also called the 'force of
#'  infection' (\eqn{\lambda}). This is different from the effective
#' transmission rate (\eqn{\beta}).
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i}.
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group \eqn{i} and risk group \eqn{k}.
#' @param p_susceptibility A matrix giving the probability that an individual
#' in demography group \eqn{i} is in risk (or susceptibility) group \eqn{k}.
#' Each row represents the overall distribution of individuals in demographic
#' group \eqn{i} across risk groups, and each row *must sum to 1.0*.
#' @param infectious_period Duration of the infectious period in days.
#' Default value is 1.8 days.
#' @export
#' @return Returns the \eqn{R_0} for the infection in the population.
#' @examples
#' # Get example dataset and prepare contact matrix and demography
#' data(polymod_uk)
#' contact_matrix <- polymod_uk$contact_matrix
#' demography_vector <- polymod_uk$demography_vector
#'
#' # define lambda
#' lambda <- 0.3
#'
#' # define infectious period of 5 days
#' infectious_period <- 5
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
#' lambda_to_r0(
#'   lambda, contact_matrix, demography_vector,
#'   susceptibility, p_susceptibility,
#'   infectious_period
#' )
#'
lambda_to_r0 <- function(lambda, contact_matrix, demography_vector,
                         susceptibility, p_susceptibility,
                         infectious_period = 1.8) {
  # input checking
  checkmate::assert_number(lambda, lower = 0, finite = TRUE)
  checkmate::assert_number(infectious_period, lower = 0, finite = TRUE)
  checkmate::assert_matrix(
    contact_matrix,
    mode = "numeric", all.missing = FALSE, any.missing = FALSE,
    min.rows = 1, min.cols = 1
  )
  checkmate::assert_numeric(
    demography_vector,
    lower = 0, finite = TRUE,
    any.missing = FALSE, all.missing = FALSE, min.len = 1
  )
  stopifnot(
    "Error: contact matrix must have as many rows as demography groups" =
      (nrow(contact_matrix) == length(demography_vector))
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
    Y = contact_matrix
  )

  # multiply rows of expanded contact matrix by susceptibility proportions
  # and multiply columns by the expanded demography vector
  contact_matrix_spread <- t(
    t(contact_matrix_spread * as.vector(susceptibility)) *
      demography_vector_spread
  )

  # calculate the largest real eigenvalue of the new matrix
  # return the product of the dominant eigenvalue, lambda, and the infectious
  # period as R0
  max(Re(eigen(contact_matrix_spread, only.values = TRUE)$values)) * lambda * infectious_period
}

#' @title Calculate transmission rate (\eqn{\lambda}) from \eqn{R_0}
#'
#' @description Uses the R0 (\eqn{R0}), contact matrix (\eqn{C}),
#' population (\eqn{N}), and infectious period (\eqn{\gamma})
#' to calculate the transmission rate using the following equation.
#' \deqn{\lambda = R_0 / ({Max}(EV(C)) \gamma)}
#' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the social contacts matrix scaled by the number of
#' individuals in each demographic and susceptibility group in the population.
#'
#' @param r0 The basic reproductive number \eqn{R_0} of the infection.
#' @inheritParams lambda_to_r0
#' @export
#' @return Returns the transmission rate of the infection, also called the
#' 'force of infection' (\eqn{\lambda}). This is different from the effective
#' transmission rate (\eqn{\beta}).
#' @examples
#' # Get example dataset and prepare contact matrix and demography
#' data(polymod_uk)
#' contact_matrix <- polymod_uk$contact_matrix
#' demography_vector <- polymod_uk$demography_vector
#'
#' # define R0 similar to pandemic influenza
#' r0 <- 1.5
#' # define infectious period of 5 days
#' infectious_period <- 5
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
#' r0_to_lambda(
#'   r0, contact_matrix, demography_vector,
#'   susceptibility, p_susceptibility,
#'   infectious_period
#' )
#'
r0_to_lambda <- function(r0, contact_matrix, demography_vector,
                         susceptibility, p_susceptibility,
                         infectious_period = 1.8) {
  # input checking
  checkmate::assert_number(r0, lower = 0, finite = TRUE)
  checkmate::assert_number(infectious_period, lower = 0, finite = TRUE)
  checkmate::assert_matrix(
    contact_matrix,
    mode = "numeric", all.missing = FALSE, any.missing = FALSE,
    min.rows = 1, min.cols = 1
  )
  checkmate::assert_numeric(
    demography_vector,
    lower = 0, finite = TRUE,
    any.missing = FALSE, all.missing = FALSE, min.len = 1
  )
  stopifnot(
    "Error: contact matrix must have as many rows as demography groups" =
      (nrow(contact_matrix) == length(demography_vector))
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
    Y = contact_matrix
  )

  # divide the contact matrix by the proportion of susceptibles,
  # and multiply by demography
  contact_matrix_spread <- t(
    t(contact_matrix_spread / as.vector(susceptibility)) *
      demography_vector_spread
  )

  # calculate the largest real eigenvalue of the new matrix
  # return lambda as the R0 divided by the (dominant eigenvalue times the
  # infectious period)
  r0 / (max(Re(eigen(contact_matrix_spread, only.values = TRUE)$values)) * infectious_period)
}
