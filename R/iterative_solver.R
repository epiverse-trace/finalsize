#' Iterative solver implemented in R
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
#' @return A vector of final sizes, of the size (N demography groups *
#' N risk groups).
solve_final_size_iterative <- function(contact_matrix,
                                       demography_vector,
                                       p_susceptibility,
                                       susceptibility,
                                       iterations = 1000,
                                       tolerance = 1e-6,
                                       step_rate = 1.9,
                                       adapt_step = TRUE) {
  # count demography groups
  n_dim <- length(demography_vector)

  # make vector of zeros
  zeros <- rep(0.0, n_dim)
  # make vector of initial final size guesses = 0.5
  epi_final_size <- rep(0.5, n_dim)

  # replicate contact matrix
  contact_matrix_ <- contact_matrix
  # set contact_matrix values to zero if there are no contacts among
  # demography groups, or if demography groups are empty
  i_here <- demography_vector == 0 | susceptibility == 0 |
    rowSums(contact_matrix) == 0
  zeros[i_here] <- 1.0
  epi_final_size[i_here] <- 0.0

  # matrix filled by columns
  contact_matrix_ <- matrix(
    as.vector(contact_matrix_) *
      (susceptibility %x% demography_vector), # note Kronecker product
    nrow = nrow(contact_matrix_),
    ncol = ncol(contact_matrix_)
  )

  contact_matrix_[i_here, zeros == 1] <- 0.0

  # make a vector to hold final size, empty numeric of size n_dim
  epi_final_size_return <- numeric(n_dim)

  # define functions to minimise error in final size estimate
  fn_f <- function(epi_final_size_, epi_final_size_return_) {
    # contact_matrix_ captured from environment
    s <- as.vector(contact_matrix_ %*% (-epi_final_size_))

    epi_final_size_return_ <- ifelse(zeros == 1.0, 0.0, 1.0)
    epi_final_size_return_ <- epi_final_size_return_ - (p_susceptibility *
      exp(susceptibility * s))

    epi_final_size_return_
  }
  # define initial current error
  current_error <- step_rate * n_dim
  # run over fn_f over epi_final_size (intial guess)
  # and epi_final_size_return (current estimate?)
  for (i in seq(iterations)) {
    epi_final_size_return <- fn_f(epi_final_size, epi_final_size_return)
    # get current error
    dpi <- epi_final_size - epi_final_size_return
    error <- sum(abs(dpi))
    # break loop if error below tolerance
    if (error < tolerance) {
      epi_final_size <- epi_final_size - dpi
      break
    }
    # adapt the step size based on the change in error
    change <- current_error - error
    if (change > 0.0) {
      epi_final_size <- epi_final_size - step_rate * dpi
      if (adapt_step) {
        step_rate <- step_rate * 1.1
      }
    } else {
      epi_final_size <- epi_final_size - dpi
      if (adapt_step) {
        step_rate <- max(0.9 * step_rate, 1.0)
      }
    }
    current_error <- error
  }

  # adjust numerical errors
  # relies on TRUE equivalent to 1
  epi_final_size <- ifelse(zeros, 0.0, epi_final_size)

  # return what
  epi_final_size
}
