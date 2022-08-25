#' Calculate final epidemic size
#'
#' This function calculates final epidemic size using SIR model for
#' heterogeneously mixing population
#' @param r0  Basic reproduction number. Default is 2.
#' @param contact_matrix  Social contact matrix. Entry mm_ij gives average
#' number of contacts in group i reported by participants in group j.
#' @param demography  Demography vector. Entry pp_i gives proportion of
#' total population in group i (model will normalise if needed)
#' @param susceptibility  Proportion of each group susceptible. Null assumption is
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
final_size_cpp <- function(r0 = 2, contact_matrix, demography,
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
  
  tryCatch(
    expr = solve_final_size_internal(
      contact_matrix = contact_matrix,
      demography = demography, 
      p_susceptibility = as.matrix(susceptibility),
      susceptibility = as.matrix(susceptibility)
    ),
    error = function(e) {
      print(e)
      message("final size cpp failed as expected")
    }
  )
}
