#' @title Final size of epidemic by susceptibility groups
#'
#' @description `final_size_grps` allows the calculation of epidemic final sizes
#' in a population with heterogeneous mixing, and with heterogeneous
#' susceptibility to infection.
#'
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
#' @param solver Which solver to use. Options are "iterative" or "newton", for
#' the iterative solver, or the Newton solver, respectively. Special conditions
#' apply when using the Newton solver, see the `control` argument.
#' @param control A list of named solver options, see *Details*.
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
#' psusc <- t(apply(psusc, 1, \(x) x / sum(x)))
#' susc <- matrix(
#'   data = 1, nrow = n_demo_grps, ncol = n_risk_grps
#' )
#'
#' control <- list(
#'   iterations = 1000,
#'   tolerance = 1e-6,
#'   step_rate = 1.9,
#'   adapt_step = TRUE
#' )
#'
#' final_size_grps(
#'   contact_matrix = r0 * c_matrix,
#'   demography_vector = d_vector,
#'   p_susceptibility = psusc,
#'   susceptibility = susc,
#'   solver = "newton"
#' )
#'
final_size_grps <- function(contact_matrix,
                            demography_vector,
                            p_susceptibility,
                            susceptibility,
                            solver = c("iterative", "newton"),
                            control = list()) {
  # check arguments input
  stopifnot(
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

  f_solver <- switch(solver,
    newton = solve_final_size_newton,
    iterative = solve_final_size_iterative
  )

  # prepare default list of control params
  con <- list(
    iterations = 1000,
    tolerance = 1e-6,
    step_rate = 1.9,
    adapt_step = TRUE
  )
  # assign user specified values
  con[(names(control))] <- control
  extra_names <- names(control)[!names(control) %in% names(con)]
  if (length(extra_names)) {
    warning("Unknown names in control: ", paste(extra_names, collapse = ", "))
  }

  # prune control list if using newton solver
  if (solver == "newton") con <- con[c("iterations", "tolerance")]

  # prepare epi spread object
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  epi_final_size <- numeric()
  # solve for final size using solver control options
  epi_final_size <- do.call(
    f_solver,
    c(
      list(
        contact_matrix = epi_spread_data[["contact_matrix"]],
        demography_vector = epi_spread_data[["demography_vector"]],
        susceptibility = epi_spread_data[["susceptibility"]]
      ),
      con
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
