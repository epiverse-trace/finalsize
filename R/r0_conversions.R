#' @title Converts between epidemiological parameters related to \eqn{R_0}
#'
#' @name r0_conversions
#' @rdname r0_conversions
#'
#' @description Converts between \eqn{R_0} and the transmission rate
#' \eqn{\lambda}, or calculates
#' the effective reproduction number \eqn{R_\text{eff}} for a population,
#' while accounting for population characteristics including demographic
#' heterogeneity in social contacts, heterogeneity in the demographic
#' distribution, and heterogeneity in susceptibility to infection.
#'
#' @description Uses the R0 (\eqn{R_0}), contact matrix (\eqn{C}),
#' population (\eqn{N}), and infectious period (\eqn{\gamma})
#' to calculate the transmission rate using the following equation.
#' \deqn{\lambda = R_0 / ({Max}(EV(C)) \gamma)}
#' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the social contacts matrix scaled by the number of
#' individuals in each demographic and susceptibility group in the population.
#'
#' @param r0 The basic reproductive number \eqn{R_0} of the infection.
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
#' @param contact_scaling For `r_eff()`, either a single number or a numeric
#' vector of the same length as `demography_vector`, giving the proportional
#' scaling of contacts of each demographic group. Values must be in the range
#' \eqn{[0, 1]}. Defaults to 1.0 for no scaling.
#'
#' @details
#' Given the transmission rate (\eqn{\lambda}),
#' social contacts matrix (\eqn{c}), demography (\eqn{N}), the
#' distribution \eqn{P} of each demographic group \eqn{i} into
#' susceptibility groups \eqn{S}, and the infectious period (\eqn{\gamma})
#'
#' - `r_eff()` calculates the effective reproductive number;
#'
#' - `lamda_to_r0()` calculates the \eqn{R_0} from the transmission rate as
#' \deqn{R_0 = {Max}(EV(C)) \times \lambda \gamma}
#'
#' - `r0_to_lambda()` calculates the transmission rate as
#' \deqn{\lambda = R_0 / ({Max}(EV(C)) \gamma)}
#' Note that this is also called the 'force of infection' and is different from
#' the effective transmission rate often denoted \eqn{\beta}.
#'
#' Here, \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is
#' calculated from the social contacts matrix scaled by the number of
#' individuals in each demographic and susceptibility group in the population.
#'
#' @export
#' @return Returns a single number for the calculated value.
#' @examples
#' #### Prepare data ####
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
#' # In this example, risk varies across groups
#' susceptibility <- matrix(
#'   data = c(0.5, 0.7, 1.0), nrow = n_demo_grps, ncol = n_risk_grps
#' )
#'
#' # risk does not vary within groups
#' p_susceptibility <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#' # p_susceptibility rows must sum to 1.0
#' p_susceptibility <- p_susceptibility / rowSums(p_susceptibility)
#'
#' #### Effective R ####
#' r0 <- 2.0
#' r_eff(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility
#' )
#'
#' # With a 5% reduction in contacts
#' r_eff(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility,
#'   contact_scaling = 0.95
#' )
#'
#' #### Transmission rate to R0 ####
#' lambda_to_r0(
#'   lambda, contact_matrix, demography_vector,
#'   susceptibility, p_susceptibility,
#'   infectious_period
#' )
#'
#' #### R0 to Transmission rate ####
#' r0 <- 1.5
#' r0_to_lambda(
#'   r0, contact_matrix, demography_vector,
#'   susceptibility, p_susceptibility,
#'   infectious_period
#' )
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
  max(Re(eigen(
    contact_matrix_spread,
    only.values = TRUE
  )$values)) * lambda * infectious_period
}

#' @title Calculate transmission rate (\eqn{\lambda}) from \eqn{R_0}
#'
#' @name r0_conversions
#' @rdname r0_conversions
#'
#' @export
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
  r0 / (max(Re(eigen(
    contact_matrix_spread,
    only.values = TRUE
  )$values)) * infectious_period)
}

#' @title Calculate \eqn{R_\text{eff}} in a heterogeneous population
#'
#' @name r0_conversions
#' @rdname r0_conversions
#'
#' @export
r_eff <- function(r0,
                  contact_matrix,
                  demography_vector,
                  susceptibility,
                  p_susceptibility,
                  contact_scaling = 1.0) {
  # check arguments input
  checkmate::assert_number(r0, lower = 0, finite = TRUE)
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
    "Error: contact_scaling must be a number or length of demography vector" =
      is.numeric(contact_scaling) && (length(contact_scaling) == 1 ||
        length(contact_scaling) == length(demography_vector)),
    "Error: contact_scaling must be in the range 0.0 -- 1.0" =
      all(contact_scaling >= 0.0 & contact_scaling <= 1.0),
    "Error: contact matrix must have a maximum real eigenvalue of 1.0" =
      (
        abs(max(Re(eigen(
          contact_matrix * demography_vector,
          only.values = TRUE
        )$values) - 1.0)) < 1e-6
      )
  )

  # count risk groups
  n_susc_groups <- ncol(p_susceptibility)
  # make lps, a vector of all p_susc values
  lps <- as.vector(p_susceptibility)

  # replicate the demography vector and multiply by p_susceptibility
  demography_vector_spread <- rep(demography_vector, n_susc_groups)
  demography_vector_spread <- demography_vector_spread * lps

  # scale both incoming and outgoing contacts, but prevent squared scaling
  n_demo_grps <- length(demography_vector)
  contact_scaling_incoming <- matrix(contact_scaling, n_demo_grps, n_demo_grps)
  contact_scaling_outgoing <- t(contact_scaling_incoming)
  diag(contact_scaling_outgoing) <- 1.0

  contact_matrix <- contact_matrix *
    contact_scaling_incoming * contact_scaling_outgoing

  # replicate contact matrix
  contact_matrix_spread <- kronecker(
    X = matrix(1, nrow = n_susc_groups, ncol = n_susc_groups),
    Y = contact_matrix * r0
  )

  # process contact matrix in solver specific ways
  contact_matrix_spread <- t(
    t(contact_matrix_spread * as.vector(susceptibility)) *
      demography_vector_spread
  )

  # get effective R
  r_eff_ <- max(Re(eigen(contact_matrix_spread, only.values = TRUE)$values))

  # return r_eff_
  r_eff_
}
