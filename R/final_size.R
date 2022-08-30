#' Calculate final epidemic size
#'
#' This function calculates final epidemic size using SIR model for
#' heterogeneously mixing population
#' @param r0  Basic reproduction number. Default is 2.
#' @param contact_matrix  Social contact matrix. Entry mm_ij gives average
#' number of contacts in group i reported by participants in group j.
#' @param demography  Demography vector. Entry pp_i gives proportion of
#' total population in group i (model will normalise if needed)
#' @param susceptibility  Proportion of each group susceptible.
#' Null assumption is fully susceptible.
#' @keywords epidemic model
#' @export
#'
#' @examples
#' # basic usage of finalsize
#' \dontrun{
#' final_size()
#' }
#'
final_size <- function(r0 = 2, contact_matrix, demography,
                       susceptibility = NULL) {

  # Check inputs
  if (is.null(susceptibility)) {
    susceptibility <- demography * 0 + 1
  } # Assume fully susceptible if no entry
  if (nrow(contact_matrix) != length(demography)) {
    stop("demography vector needs to be same size as contact matrix")
  }
  if (length(demography) != length(susceptibility)) {
    stop("demography vector needs to be same size as susceptibility vector")
  }
  pp0 <- as.numeric(demography / sum(demography))

  # Scale next generation matrix so max eigenvalue=r0
  mm0 <- as.matrix(contact_matrix)
  mm0 <- r0 * mm0 / max(Re(eigen(mm0)$values))

  # Define transmission matrix A = mm_ij*pp_j/pp_i
  beta1 <- (mm0 * susceptibility) / pp0
  beta2 <- t(t(beta1) * pp0)

  # Newton method for solving final size equation: A(1-x) = -log(x)
  vsize <- length(pp0)

  f1 <- function(beta2, x) {
    beta2 %*% (1 - x) + log(x)
  }
  f2 <- function(beta2, x, size) {
    -beta2 + diag(size) / x
  }

  # Set storage vector and precision
  iterations <- 30
  iterate_output <- matrix(NA_real_, nrow = iterations, ncol = vsize)
  x0 <- 0.001 * pp0 # Set starting point

  for (ii in seq(iterations)) {
    if (ii == 1) {
      xx <- x0
      iterate_output[ii, ] <- xx
    } else {
      dx <- solve(f2(beta2, xx, vsize), -f1(beta2, xx))
      iterate_output[ii, ] <- xx + dx
      xx <- as.numeric(xx + dx)
    }
  }

  (1 - iterate_output[iterations, ])
}
