#' @title Final size of an epidemic
#'
#' @description `final_size` calculates the final size of an epidemic outbreak
#' in a population with heterogeneous mixing, and with heterogeneous
#' susceptibility to infection such as that conferred by an immunisation
#' programme.
#'
#' @details
#' # Solver options
#'
#' The `control` argument accepts a list of solver options, with the iterative
#' solver taking two extra arguments than the Newton solver.
#'
#' ## Common options
#'
#' 1. `iterations`: The number of iterations over which to solve for the final
#' size, unless the error is below the solver tolerance.
#' 2. `tolerance`: The solver tolerance, set to `1e-6` by default; solving for
#' final size ends when the error drops below this tolerance.
#'
#' ## Iterative solver options
#' 1. `step_rate`: The solver step rate. Defaults to 1.9 as a value found to
#' work well.
#' 2. `adapt_step`: Boolean, whether the solver step rate should be changed
#' based on the solver error. Defaults to TRUE.
#'
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i} (model will normalise
#' if needed)
#' @param p_susceptibility A matrix giving the probability that an individual
#' in demography group \eqn{i} is in risk (or susceptibility) group \eqn{j}.
#' Each row represents the overall distribution of individuals in demographic
#' group $i$ across risk groups, and each row *must sum to 1.0*.
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group \eqn{i} and risk group \eqn{j}.
#' @param solver Which solver to use. Options are "iterative" (default) or
#' "newton", for the iterative solver, or the Newton solver, respectively.
#' Special conditions apply when using the Newton solver, see the `control`
#' argument.
#' @param control A list of named solver options, see *Solver options*.
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
#'   age.limits = c(0, 20, 40),
#'   split = TRUE # scaling by demography
#' )
#' c_matrix <- t(contact_data$matrix)
#' d_vector <- contact_data$participants$proportion
#' n_demo_grps <- length(d_vector)
#' n_risk_grps <- 3
#' # prepare p_susceptibility and susceptibility
#' psusc <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#' psusc <- psusc / rowSums(psusc)
#' susc <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#'
#' # using default settings - this selects the iterative solver
#' final_size(
#'   contact_matrix = r0 * c_matrix,
#'   demography_vector = d_vector,
#'   p_susceptibility = psusc,
#'   susceptibility = susc
#' )
#'
#' # using manually specified solver settings
#' control <- list(
#'   iterations = 1000,
#'   tolerance = 1e-6,
#'   step_rate = 1.9,
#'   adapt_step = TRUE
#' )
#'
#' final_size(
#'   contact_matrix = r0 * c_matrix,
#'   demography_vector = d_vector,
#'   p_susceptibility = psusc,
#'   susceptibility = susc,
#'   solver = "iterative",
#'   control = control
#' )
#'
#' # manual settings for the newton solver
#' control <- list(
#'   iterations = 50000,
#'   tolerance = 1e-12
#' )
#'
#' final_size(
#'   contact_matrix = r0 * c_matrix,
#'   demography_vector = d_vector,
#'   p_susceptibility = psusc,
#'   susceptibility = susc,
#'   solver = "newton",
#'   control = control
#' )
#'
final_size <- function(contact_matrix,
                       demography_vector,
                       p_susceptibility,
                       susceptibility,
                       solver = c("iterative", "newton"),
                       control = list(
                         iterations = 10000,
                         tolerance = 1e-6,
                         step_rate = 1.9,
                         adapt_step = TRUE
                       )) {
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
      )
  )

  # check which solver is requested
  solver <- match.arg(arg = solver, several.ok = FALSE)

  epi_final_size <- do.call(
    final_size_,
    c(
      list(
        contact_matrix,
        demography_vector,
        p_susceptibility,
        susceptibility,
        solver
      ),
      control
    )
  )

  # unroll p_susceptibility data
  lps <- as.vector(p_susceptibility)

  # multiply demo-risk specific final sizes by corresponding pop proportions
  epi_final_size <- epi_final_size * lps

  # final sizes mapped to matrix with dims (n_demo_grp, n_risk_grps)
  epi_final_size <- matrix(
    data = epi_final_size,
    nrow = length(demography_vector),
    ncol = ncol(p_susceptibility)
  )
  # return row-wise sum, i.e., the demo-grp wise sum
  rowSums(epi_final_size)
}
