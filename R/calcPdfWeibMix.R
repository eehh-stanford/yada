
#' @title Calculate PDF given Weibull mixture proportions and parameters
#'
#' @description Given the input vector x and Weibull mixture proportions and parameters z and theta, calculate the probability density function for the input vector x.
#'
# @details
#'
#' @param x Locations at which to evaluate probability density function
#' @param z Vector of mixture proportions
#' @param theta The value of shape and scale parameters with the ordering [sh1,sc1,sh2,sc2,...]
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The PDF vector for the input vector x
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{fitWeibMixMaxLikEM,calcPdfMatrixWeibMix}
# @references

calcPdfWeibMix <- function(x, z, theta) {
  pdfMatrix <- calcPdfMatrixWeibMix(x, theta)
  pdf <- rowSums(t(t(pdfMatrix) * z))
  return(pdf)
}
