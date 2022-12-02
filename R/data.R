#' Example POLYMOD social contact data for the U.K.
#'
#' An example of social contact and demography data for use with `finalsize`,
#' accessed from the POLYMOD social contacts dataset using the `socialmixr`
#' package. Data are for the United Kingdom, and age limits are set at 0, 20,
#' and 40 years, with \code{symmetric = TRUE}. Code to get these data is
#' given in `data-raw/polymod_uk.R`.
#'
#' @format ## `polymod_uk`
#' A list with two named elements:
#' \describe{
#'   \item{contact_matrix}{A contact matrix with mean contacts between age
#'      groups. This matrix is scaled by its largest real eigenvalue, and  each
#'      row is scaled by the corresponding element in the `demography_vector`.}
#'   \item{demography_vector}{A vector with the number of individuals in each of
#'      three age groups: 0 -- 20, 20 -- 40, 40+.}
#' }
#' @source \doi{10.1371/journal.pmed.0050074}; obtained using
#' \code{socialmixr::polymod}. See further methods in `data-raw/polymod_uk.R`.
"polymod_uk"
