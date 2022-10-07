#' Newton solver implemented in R
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param susceptibility A matrix giving the susceptibility of individuals in
#' demographic group $i$ and risk group $j$.
#' @param iterations Number of solver iterations. Defaults to 1,000.
#' @param tolerance Solver error tolerance.
#'
#' @return A vector final sizes, of the size (N demography groups *
#' N risk groups).
solve_final_size_newton <- function(contact_matrix,
                                    demography_vector,
                                    susceptibility,
                                    iterations = 1000,
                                    tolerance = 1e-6) {
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
  contact_matrix_ <- contact_matrix * demography_vector %o% susceptibility

  contact_matrix_[i_here, i_here] <- 0.0

  # make a copy of the initial guesses for final size
  cache_v <- epi_final_size

  # define functions for matrix operations
  f1 <- function(x) {
    # contact_matrix_ captured from environment
    (contact_matrix_ %*% (1.0 - x)) + log(x)
  }

  f2 <- function(x) {
    # contact_matrix_ captured from environment
    -contact_matrix_ + diag(1.0 / x)
  }

  # cache_m is the contact matrix. must not be modified!
  dx_f <- function(x, cache, cache_m) {
    cache_m <- f2(x)
    cache <- -f1(x)
    # partial pivoting LU decomposition
    cache_m_pivlu <- Matrix::lu(cache_m)
    cache_m_pivlu <- Matrix::expand(cache_m_pivlu)
    cache_m_pivlu <- cache_m_pivlu$L %*% cache_m_pivlu$U
    # return solution
    solve(cache_m_pivlu, cache)
  }

  # prepare an initial value from where to
  x_ <- rep(1e-6, n_dim)

  for (i in seq(iterations)) {
    cache_v <- as.vector(dx_f(x_, cache_v, contact_matrix_))
    error_ <- sum(abs(cache_v))
    x_ <- x_ + cache_v

    if (error_ < tolerance) {
      break
    }
  }

  epi_final_size <- 1 - x_ # returns a vector
  epi_final_size
}
