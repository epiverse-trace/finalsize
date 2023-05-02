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
#' error drops below this tolerance. Defaults to set `1e-6`. Larger tolerance
#' values are likely to lead to inaccurate final size estimates.
#'
#' ## Iterative solver options
#' 1. `step_rate`: The solver step rate. Defaults to 1.9 as a value found to
#' work well.
#' 2. `adapt_step`: Boolean, whether the solver step rate should be changed
#' based on the solver error. Defaults to TRUE.
#'
#' @param r0 The basic reproductive number \eqn{R_0} of the disease.
#' @param contact_matrix Social contact matrix. Entry \eqn{m_{ij}} gives average
#' number of contacts in group \eqn{i} reported by participants in group \eqn{j}
#' @param demography_vector Demography vector. Entry \eqn{v_{i}} gives
#' proportion of total population in group \eqn{i} (model will normalise
#' if needed).
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group \eqn{i} and risk group \eqn{k}.
#' @param p_susceptibility A matrix giving the probability that an individual
#' in demography group \eqn{i} is in risk (or susceptibility) group \eqn{k}.
#' Each row represents the overall distribution of individuals in demographic
#' group \eqn{i} across risk groups, and each row *must sum to 1.0*.
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
#' # using default arguments for `solver` and `control`
#' final_size(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility
#' )
#'
#' # using manually specified solver settings for the iterative solver
#' control <- list(
#'   iterations = 100,
#'   tolerance = 1e-3,
#'   step_rate = 1.9,
#'   adapt_step = TRUE
#' )
#'
#' final_size(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility,
#'   solver = "iterative",
#'   control = control
#' )
#'
#' # manual settings for the newton solver
#' control <- list(
#'   iterations = 100,
#'   tolerance = 1e-3
#' )
#'
#' final_size(
#'   r0 = r0,
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   susceptibility = susceptibility,
#'   p_susceptibility = p_susceptibility,
#'   solver = "newton",
#'   control = control
#' )
final_size <- function(r0,
                       contact_matrix,
                       demography_vector,
                       susceptibility,
                       p_susceptibility,
                       solver = c("iterative", "newton"),
                       control = list()) {
  # check arguments input
  checkmate::assert_number(r0, lower = 0, finite = TRUE)
  stopifnot(
    "Error: contact matrix must be a matrix" =
      (is.matrix(contact_matrix)),
    "Error: demography vector must be a numeric vector" =
      (is.vector(demography_vector, mode = "numeric"))
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
      ),
    "Error: control list names can only be: 'iterations', 'tolerance',
    'step_rate', 'adapt_step'" =
      (
        all(names(control) %in% c(
          "iterations", "tolerance",
          "step_rate", "adapt_step"
        )) || (length(control) == 0)
      )
  )

  # prepare default control list
  con <- list(
    iterations = 10000,
    tolerance = 1e-6,
    step_rate = 1.9,
    adapt_step = TRUE
  )
  # pass user solver options to default control list
  con[names(control)] <- control

  # # check which solver is requested
  solver <- match.arg(arg = solver, several.ok = FALSE)

  # prepare the population data for the solver
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
    Y = r0 * contact_matrix
  )

  # get group wise final sizes
  epi_final_size <- .final_size(
    c(
      list(
        contact_matrix = contact_matrix_spread,
        demography_vector = demography_vector_spread,
        susceptibility = as.vector(susceptibility),
        solver = solver
      ),
      con
    )
  ) # using internal Rcpp function

  # check for names of age groups and susceptibility groups
  names_demography <- names(demography_vector)
  if (is.null(names_demography)) {
    # check if contact matrix has column names
    if (!is.null(colnames(contact_matrix))) {
      names_demography <- colnames(contact_matrix)
    } else if (!is.null(rownames(contact_matrix))) {
      names_demography <- rownames(contact_matrix)
    } else {
      names_demography <- sprintf(
        "demo_grp_%i", seq_along(demography_vector)
      )
    }
  }

  # check for susceptibility group names
  names_susceptibility <- colnames(susceptibility)
  if (is.null(names_susceptibility)) {
    names_susceptibility <- sprintf(
      "susc_grp_%i", seq_len(ncol(susceptibility))
    )
  }

  epi_final_size <- data.frame(
    demo_grp = rep(
      names_demography,
      times = ncol(susceptibility)
    ),
    susc_grp = rep(
      names_susceptibility,
      each = nrow(susceptibility)
    ),
    susceptibility = as.vector(susceptibility),
    p_infected = epi_final_size
  )
  epi_final_size
}
