#' Calculate final epidemic size
#'
#' @param contact_matrix 
#' @param demography 
#' @param susceptibility 
#' @param p_susceptibility 
#' @param solver_adapt_step 
#' @param solver_tolerance 
#' @return A matrix of final proportions of age and susceptibility classes
#' infected. WIP.
#' 
#' @export
#' 
final_size_cpp <- function(contact_matrix, demography,
                           susceptibility, p_susceptibility,
                           solver_adapt_step = TRUE,
                           solver_tolerance = 1e-6) {
  # need checks on arguments here
  # scale demography proportions to 1
  demography = demography / sum(demography)
  
  contact_matrix <- t(contact_matrix$matrix)
  # scale contact matrix cols by prop of dem
  # for(i in seq(ncol(contact_matrix))) contact_matrix[, i] = contact_matrix[, i] / demography[i]
  
  # scale contact matrix by largest real eigenvalue of the matrix
  contact_matrix = contact_matrix / (max(Re(eigen(contact_matrix)$values)))
  
  tryCatch(
    expr = solve_final_size_internal(
      contact_matrix = contact_matrix,
      demography = demography,
      p_susceptibility = p_susceptibility,
      susceptibility = susceptibility
    ),
    error = function(e) {
      print(e)
      message("final size cpp failed as expected")
    }
  )
}
