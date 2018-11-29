#' @title Calculate the negative log-likelihood given Weibull mixture proportions and parameters
#'
#' @description Given the input vector x and Weibull mixture proportions and parameters z and theta, calculate the negative log-likelihood.
#'
# @details
#'
#' @param theta The value of shape and scale parameters with the ordering [sh1,sc1,sh2,sc2,...]
#' @param x Locations at which to evaluate negative log-likelihood
#' @param z Vector of mixture proportions
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The negative log-likelihood (summed)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{fitWeibMixMaxLikEM,calcPdfWeibMix}
# @references

calcNegLogLikWeibMix <- function(theta, x, z) {
  f <- calcPdfWeibMix(x, z, theta)
  nll <- -sum(log(f))
  return(nll)
}
