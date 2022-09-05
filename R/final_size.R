#' Calculate final epidemic size
#'
#' This function calculates final epidemic size using SIR model for
#' heterogeneously mixing population
#' @param r0  Basic reproduction number. Default is 2
#' @param contact_matrix  Social contact matrix. Entry $mm_{ij}$ gives average
#' number of contacts in group $i$ reported by participants in group j
#' @param demography_vector  Demography vector. Entry $pp_{i}$ gives proportion
#' of total population in group $i$ (model will normalise if needed)
#' @param prop_suscep  Proportion of each group susceptible. Null assumption is
#' fully susceptible
#' @keywords epidemic model
#' @export
#' @examples
#' library(socialmixr)
#' data(polymod)
#' contact_data <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0,20,40))
#' 
#' c_matrix <-t(contact_data$matrix) # Define contact matrix (entry ij is contacts in group i reported by group j)
#' d_vector <- contact_data$participants$proportion # Define proportion in each age group
#' p_suscep <- c(1,0.5,0.5) # Define proportion of age group that is susceptible to infection
#' 
#' # Run final size model
#' final_size(r0, contact_matrix = c_matrix, demography_vector = d_vector, prop_suscep = p_suscep)
#'
final_size <- function(r0 = 2, contact_matrix, demography_vector,
                       prop_suscep = NULL) {

  # Check inputs
  if (is.null(prop_suscep)) {
    prop_suscep <- rep(1, length(demography_vector))
  } # Assume fully susceptible if no entry
  
  checkmate::assert_count(r0)
  checkmate::assert_vector(demography_vector)
  checkmate::assert_matrix(contact_matrix)
  checkmate::assert_atomic(prop_initial_infected)
  
  assertthat::assert_that(
    nrow(contact_matrix) == length(demography_vector),
    msg = "demography vector needs to be same size as contact matrix"
  )
  assertthat::assert_that(
    length(demography_vector) == length(prop_suscep),
    msg = "demography vector needs to be same size as susceptibility vector"
  )
  
  }
  pp0 <- as.numeric(demography_vector / sum(demography_vector))

  # Scale next generation matrix so max eigenvalue=r0
  mm0 <- as.matrix(contact_matrix)
  mm0 <- r0 * mm0 / max(Re(eigen(mm0)$values))

  # Define transmission matrix A = mm_ij*pp_j/pp_i
  beta1 <- (mm0 * prop_suscep) / pp0
  beta2 <- t(t(beta1) * pp0)

  # Newton method for solving final size equation: A(1-x) = -log(x)

  # Define functions f and f'
  vsize <- length(pp0)
  f1 <- function(beta2, x) {
    beta2 %*% (1 - x) + log(x)
  }
  f2 <- function(beta2, x, size) {
    -beta2 + diag(size) / x
  }

  # Set storage vector and precision
  iterations <- 30
  iterate_output <- matrix(NA, nrow = iterations, ncol = vsize)
  x0 <- 0.001 * pp0 # Set starting point

  for (ii in 1:iterations) {
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
