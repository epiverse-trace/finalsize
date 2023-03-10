#' @title Calculate \eqn{R_0} from transmission rate (\eqn{\lambda})
#'
#' @description Uses the transmission rate (\eqn{\lambda}), contact matrix
#' (\eqn{c}), population (\eqn{N}), and infectious period (\eqn{\gamma}) to
#' calculate the R0 using the following equation.
#' \deqn{\lambda max(EV(C)) \gamma}
#' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the contact matrix and the population
#' (\eqn{C[i,j] = c[i,j] N[j]}).
#'
#' @param lambda The transmission rate of the disease, also called the 'force of
#'  infection' (\eqn{\lambda}). This is different from the effective
#' transmission rate (\eqn{\beta}).
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i}.
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
#'
#' lambda_to_r0(
#'   lambda, contact_matrix, demography_vector, infectious_period
#' )
#'
lambda_to_r0 <- function(lambda, contact_matrix, demography_vector,
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

  # multiply each column of the contact matrix with the corresponding
  # demography vector element
  cm <- t(t(contact_matrix) * demography_vector)

  # calculate the largest real eigenvalue of the new matrix
  # return the product of the dominant eigenvalue, lambda, and the infectious
  # period as R0
  max(Re(eigen(cm)$values)) * lambda * infectious_period
}

#' @title Calculate \eqn{R_0} from transmission rate (\eqn{\lambda})
#'
#' @description Uses the transmission rate (\eqn{\lambda}), contact matrix
#' (\eqn{c}), population (\eqn{N}), and infectious period (\eqn{\gamma}) to
#' calculate the R0 using the following equation.
#' \deqn{\lambda max(EV(C)) \gamma}
#' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the contact matrix and the population
#' (\eqn{C[i,j] = c[i,j] N[j]}).
#'
#' @param r0 The basic reproductive number \eqn{R_0} of the infection.
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i}.
#' @param infectious_period Duration of the infectious period in days.
#' Default value is 1.8 days.
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
#'
#' # define infectious period of 5 days
#' infectious_period <- 5
#'
#' r0_to_lambda(
#'   r0, contact_matrix, demography_vector, infectious_period
#' )
#'
r0_to_lambda <- function(r0, contact_matrix, demography_vector,
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

  # multiply each column of the contact matrix with the corresponding
  # demography vector element
  cm <- t(t(contact_matrix) * demography_vector)

  # calculate the largest real eigenvalue of the new matrix
  # return lambda as the R0 divided by the (dominant eigenvalue times the
  # infectious period)
  r0 / (max(Re(eigen(cm)$values)) * infectious_period)
}
