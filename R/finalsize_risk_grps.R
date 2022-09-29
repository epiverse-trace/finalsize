#' Final size of epidemic by susceptibility groups
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
#' @param iterations Number of solver iterations. Defaults to 1,000.
#' @param tolerance Solver error tolerance.
#' @param step_rate The solver step rate. Defaults to 1.9 as a value found to
#' work well.
#' @param adapt_step Boolean, whether the solver step rate should be changed
#' based on the solver error. Defaults to TRUE.
#'
#' @keywords epidemic model
#' @return A vector of final sizes by demography group.
#' @export
#' @examples
#' library(socialmixr)
#' data(polymod)
#' r0 <- 2.0
#' contact_data <- contact_matrix(
#'   polymod,
#'   countries = "United Kingdom",
#'   age.limits = c(0, 20, 40)
#' )
#' c_matrix <- t(contact_data$matrix)
#' d_vector <- contact_data$participants$proportion
#' # Scale contact matrix to demography
#' c_matrix <- apply(
#'   c_matrix, 1, function(r) r / d_vector
#' )
#' n_demo_grps <- length(d_vector)
#' n_risk_grps <- 3
#' # prepare p_susceptibility and susceptibility
#' psusc <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#' psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
#' susc <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#' final_size_grps(
#'   contact_matrix = r0 * c_matrix,
#'   demography_vector = d_vector,
#'   p_susceptibility = psusc,
#'   susceptibility = susc
#' )
#'
final_size_grps <- function(contact_matrix,
                            demography_vector,
                            p_susceptibility,
                            susceptibility,
                            iterations = 1000,
                            tolerance = 1e-6,
                            step_rate = 1.9,
                            adapt_step = TRUE) {
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
  # p_susceptibility must have rowwise sums of 1.0 or nearly 1.0
  stopifnot(
    "Error: p_susceptibility rows must sum to 1.0" =
      (
        all(abs(rowSums(p_susceptibility) - 1) < 1e-6)
      )
  )

  # prepare epi spread object
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  # solve for final size
  pi <- solve_final_size_iterative(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility_"]],
    susceptibility = epi_spread_data[["susceptibility"]],
    iterations = iterations,
    tolerance = tolerance
  )

  # unroll p_susceptibility data
  lps <- as.vector(p_susceptibility)

  # multiply demo-risk specific final sizes by corresponding pop proportions
  pi <- pi * lps

  # final sizes mapped to matrix with dims (n_demo_grp, n_risk_grps)
  pi <- matrix(
    data = pi,
    nrow = length(demography_vector),
    ncol = ncol(p_susceptibility)
  )
  # return row-wise sum, i.e., the demo-grp wise sum
  rowSums(pi)
}
