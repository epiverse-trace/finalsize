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
  pi <- rep(0.5, n_dim)

  # replicate contact matrix
  contact_matrix_ <- contact_matrix
  # set some values to zero if there are no contacts among
  # demography groups, or if demography groups are empty
  for (i in seq(nrow(contact_matrix_))) {
    if (demography_vector[i] == 0 || susceptibility[[i]] == 0 ||
      sum(contact_matrix[i, ]) == 0) {
      zeros[i] <- 1.0
      pi[i] <- 0.0
    }

    for (j in seq(ncol(contact_matrix))) {
      if (zeros[j] == 1) {
        contact_matrix_[i, j] <- 0.0
      } else {
        contact_matrix_[i, j] <- susceptibility[i] *
          contact_matrix_[i, j] * demography_vector[j]
      }
    }
  }
  # make a vector to hold final size, empty numeric of size n_dim
  pi_return <- numeric(n_dim)

  # define functions to minimise error in final size estimate
  fn_f <- function(pi_, pi_return_) {
    s <- contact_matrix_ %*% (-pi_) # contact_matrix_ captured from environment
    for (i in seq(nrow(contact_matrix_))) {
      if (zeros[i] == 1.0) {
        pi_return_[i] <- 0
        next
      }
      pi_return_[i] <- 1

      for (k in seq(ncol(p_susceptibility))) {
        pi_return_[i] <- pi_return_[i] - (p_susceptibility[i, k]) *
          exp(susceptibility[i, k] * s[i])
      }
    }
    pi_return_
  }
  # define initial current error
  current_error <- step_rate * n_dim
  # run over fn_f over pi (intial guess) and pi_return (current estimate?)
  for (i in seq(iterations)) {
    pi_return <- fn_f(pi, pi_return)
    # get current error
    dpi <- pi - pi_return
    error <- sum(abs(dpi))
    # break loop if error below tolerance
    if (error < tolerance) {
      pi <- pi - dpi
      break
    }
    # adapt the step size based on the change in error
    change <- current_error - error
    if (change > 0.0) {
      pi <- pi - step_rate * dpi
      if (adapt_step) {
        step_rate <- step_rate * 1.1
      }
    } else {
      pi <- pi - dpi
      if (adapt_step) {
        step_rate <- max(0.9 * step_rate, 1.0)
      }
    }
    current_error <- error
  }

  # adjust numerical errors
  for (i in seq(length(pi))) {
    if (zeros[i]) {
      pi[i] <- 0
    } # relies on FALSE equivalent to 0
  }

  # return what
  pi
}
