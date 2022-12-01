#' Example POLYMOD social contact data for the U.K.
#'
#' An example of social contact and demography data for use with `finalsize`,
#' accessed from the POLYMOD social contacts dataset using the `socialmixr`
#' package. Data are for the United Kingdom, and age limits are set at 0, 20,
#' and 40 years, with \code{symmetric = TRUE}.
#'
#' @format ## `polymod_uk`
#' A list with three named elements:
#' \describe{
#'   \item{matrix}{A contact matrix with mean contacts between age groups}
#'   \item{demography}{A data.frame with the populations and proportions of age
#'      groups in the U.K. population in 2005}
#'   \item{participants}{A data.frame with the number and proportion of
#'      participants in each age group in the POLYMOD survey}
#' }
#' @source \doi{10.1371/journal.pmed.0050074}; obtained using
#' \code{socialmixr::polymod}.
"polymod_uk"
