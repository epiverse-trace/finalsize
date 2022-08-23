#' Calculate final epidemic size
#'
#' This function calculates final epidemic size using SIR model for
#' heterogeneously mixing population
#' @param r0  Basic reproduction number. Default is 2.
#' @param contact_matrix  Social contact matrix. Entry mm_ij gives average
#' number of contacts in group i reported by participants in group j.
#' @param demography_vector  Demography vector. Entry pp_i gives proportion of
#' total population in group i (model will normalise if needed)
#' @param prop_suscep  Proportion of each group susceptible. Null assumption is
#' fully susceptible.
#' @keywords epidemic model
#' @export
#'
#' @examples
#' # basic usage of finalsize
#' \dontrun{
#' final_size()
#' }
#'
final_size <- function(r0 = 2, contact_matrix, demography_vector,
                       prop_suscep = NULL) {

  # Check inputs
  if (is.null(prop_suscep)) {
    prop_suscep <- demography_vector * 0 + 1
  } # Assume fully susceptible if no entry
  if (nrow(contact_matrix) != length(demography_vector)) {
    stop("demography vector needs to be same size as contact matrix")
  }
  if (length(demography_vector) != length(prop_suscep)) {
    stop("demography vector needs to be same size as susceptibility vector")
  }
  pp0 <- as.numeric(demography_vector / sum(demography_vector))

  # Scale next generation matrix so max eigenvalue=r0
  mm0 <- as.matrix(contact_matrix)
  mm0 <- r0 * mm0 / max(Re(eigen(mm0)$values))

  # Define transmission matrix A = mm_ij*pp_j/pp_i
  beta1 <- (mm0 * prop_suscep) / pp0
  beta2 <- t(t(beta1) * pp0)

  # Newton method for solving final size equation: A(1-x) = -log(x)
  vsize <- length(pp0)

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
