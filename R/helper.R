# Define functions f and f'
#' Title
#'
#' @param beta2 
#' @param x 
#'
f1 <- function(beta2, x) {
    beta2 %*% (1 - x) + log(x)
}

#' Title
#'
#' @param beta2 
#' @param x 
#' @param size 
#'
f2 <- function(beta2, x, size) {
    -beta2 + diag(size) / x
}