#' @title Calculate PDF matrix given Weibull mixture parameters
#'
#' @description Given the input vector x and Weibull mixture parameters theta, calculate the probability density function (PDF) matrix, which has dimensions length(x) by numMix, where numMix = length(theta)/2. The mixture proportions are not accounted for in the calculation.
#'
# @details
#'
#' @param x Locations at which to evaluate probability density function
#' @param theta The value of shape and scale parameters with the ordering [sh1,sc1,sh2,sc2,...]
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The PDF matrix with dimensions length(x) by length(theta)/2
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{fitWeibMixMaxLikEM,calcPdfWeibMix}
# @references

calcPdfMatrixWeibMix <- function(x, theta) {
  numObs <- length(x) # Number of observations
  numMix <- length(theta) / 2 # Number of mixtures
  # For x, the shape parameter, and scale parameter create a long vector of
  # length numObs*numMix for vectorized input to dweibull
  x_vect <- rep(x, numMix)
  shape_vect <- as.vector(t(matrix(rep(theta[seq(1, length(theta), by = 2)], length(x)), nrow = numMix)))
  scale_vect <- as.vector(t(matrix(rep(theta[seq(2, length(theta), by = 2)], length(x)), nrow = numMix)))
  # Calculate the probability density and transform to a matrix with dimensions numObs x numMix
  return(matrix(dweibull(x_vect, shape_vect, scale_vect), ncol = numMix))
}
