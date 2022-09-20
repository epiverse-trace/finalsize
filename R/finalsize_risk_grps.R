
#' @title Prepare contact matrix data for multiple risk groups
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param p_susceptibility WIP - what is the correct definition.
#' @param susceptibility WIP - what is the correct definition.
#'
#' @return
epi_spread <- function(contact_matrix,
                       demography_vector,
                       p_susceptibility,
                       susceptibility) {
  # replicating epi_spread from finalsize.h # --- make function later
  n_susc_groups <- ncol(p_susceptibility)
  # make p_susceptibility matrix of ones
  p_susceptibility_ <- matrix(
    1.0,
    nrow = prod(dim(p_susceptibility)), ncol = 1
  )
  # make lps, a 1 col matrix of all p_susc values
  lps <- matrix(p_susceptibility, nrow = prod(dim(p_susceptibility)), ncol = 1)
  # replicate the demography vector and multiply by p_susceptibility
  demography_vector_ <- rep(demography_vector, n_susc_groups)
  demography_vector <- demography_vector * as.vector(lps)

  # replicate contact matrix
  contact_matrix_ <- do.call("rbind", rep(list(contact_matrix), n_susc_groups))
  contact_matrix_ <- do.call(
    "cbind", rep(list(contact_matrix_), n_susc_groups)
  )

  # unroll the susceptibility matrix
  susceptibility_ <- matrix(
    data = as.vector(susceptibility),
    nrow = prod(dim(susceptibility)), ncol = 1
  )

  list(
    contact_matrix = contact_matrix_,
    demography_vector = demography_vector_,
    p_susceptibility_ = p_susceptibility_,
    susceptibility = susceptibility_
  )
}

#' Newton solver implemented in R
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param p_susceptibility WIP - what is the correct definition.
#' @param susceptibility WIP - what is the correct definition.
#' @param iterations Number of solver iterations
#' @param tolerance Solver error tolerance.
#'
#' @return
solve_final_size_newton <- function(contact_matrix,
                                    demography_vector,
                                    p_susceptibility,
                                    susceptibility,
                                    iterations = 1000
                                    tolerance = 1e-6) {
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
  # make a copy of the initial guesses for final size
  cache_v <- pi

  # define functions for matrix operations
  f1 <- function(c_m, x, cache) {
    cache <- (c_m * (1.0 - x)) + log(x)
    cache
  }

  f2 <- function(c_m, x, cache) {
    cache <- diag(1.0 / x)
    cache <- -c_m + cache
    cache
  }

  dx_f <- function(c_m, x, cache) {
    cache_m <- f2(c_m, x, cache)
    cache <- -f1(c_m, x, cache)
    cache <- solve(cache_m, cache)
    cache
  }

  # prepare inverse of the guessed final size
  x_ <- 1.0 - pi

  for (i in seq(iterations)) {
    cache_v <- dx_f(contact_matrix_, x_, cache_v)
    error_ <- sum(abs(cache_v))
    x_ <- x_ + cache_v

    if (error_ < tolerance) {
      break
    }
  }

  pi <- 1 - x_ # returns a mtrix
  pi
}

#' Final size of epidemic by susceptibility groups
#'
#' @param contact_matrix Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param p_susceptibility WIP - what is the correct definition.
#' @param susceptibility WIP - what is the correct definition.
#' @param iterations Number of solver iterations
#' @param tolerance Solver error tolerance.
#' 
#' @return A vector of final sizes by demography group.
#' @export
#' @example 
final_size_grps <- function(contact_matrix,
                                               demography_vector,
                                               p_susceptibility,
                                               susceptibility,
                                               iterations = 1000,
                                               tolerance = 1e-6) {
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

  # prepare epi spread object
  epi_spread_data <- epi_spread(
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    p_susceptibility = p_susceptibility,
    susceptibility = susceptibility
  )

  # solve for final size
  pi <- solve_final_size_newton(
    contact_matrix = epi_spread_data[["contact_matrix"]],
    demography_vector = epi_spread_data[["demography_vector"]],
    p_susceptibility = epi_spread_data[["p_susceptibility"]],
    susceptibility = epi_spread_data[["susceptibility"]],
    iterations = iterations,
    tolerance = tolerance
  )

  pi_2 <- matrix(pi,
    nrow = length(demography_vector),
    ncol = ncol(p_susceptibility)
  )

  lps <- matrix(p_susceptibility, nrow = prod(dim(p_susceptibility)), ncol = 1)
  v <- as.vector(pi_2) * lps

  pi_3 <- matrix(
    v,
    nrow = length(demography_vector), ncol = ncol(p_susceptibility)
  )

  pi_3
}
