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
#' solver taking two extra arguments than the Newton solver. This is an optional
#' argument, and default options are used within the solver functions if an
#' argument is missing. Arguments provided override the solver defaults.
#'
#' ## Common options
#'
#' 1. `iterations`: The number of iterations over which to solve for the final
#' size, unless the error is below the solver tolerance. Default = 10000.
#' 2. `tolerance`: The solver tolerance; solving for final size ends when the
#' error drops below this tolerance. Defaults to set `1e-6`.
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
#' @return A data.frame of the proportion of infected individuals, within each
#' demography group and susceptibility group combination.
#' If the demography groups and susceptibility groups are named, these
#' names are added to relevant columns. If the groups are not named, synthetic
#' names are added (e.g. `demo_grp_1`, `susc_grp_1`).
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
                       control = list()) {
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
  fn_solver <- switch(solver,
    iterative = .solve_iterative,
    newton = .solve_newton
  )

  # prepare the population data for the solver
  epi_spread_data <- .epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  # get group wise final sizes
  epi_final_size <- do.call(
    fn_solver,
    c(
      epi_spread_data,
      control
    )
  )

  # check for names of age groups and susc groups
  names_demography <- names(demography_vector)
  if (is.null(names(names_demography))) {
    # check if contact matrix has rownames
    if (!is.null(rownames(contact_matrix))) {
      names_demography <- rownames(contact_matrix)
    } else {
      names_demography <- sprintf(
        "demo_grp_%i", seq_len(length(demography_vector))
      )
    }
  }

  # check for susc group names
  names_susceptibility <- colnames(susceptibility)
  if (is.null(names(names_susceptibility))) {
    names_susceptibility <- sprintf(
      "susc_grp_%i", seq_len(ncol(susceptibility))
    )
  }

  # multiply final size by p_susceptibility
  epi_final_size <- p_susceptibility * epi_final_size
  epi_final_size <- data.frame(
    demo_grp = rep(
      names_demography,
      times = ncol(epi_final_size)
    ),
    susc_grp = rep(
      names_susceptibility,
      each = nrow(epi_final_size)
    ),
    susceptibility = as.vector(susceptibility),
    p_infected = as.vector(epi_final_size)
  )
  epi_final_size
}
