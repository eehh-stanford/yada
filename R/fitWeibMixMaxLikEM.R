#' @title Maximum likelihood fit for Weibull mixture
#'
#' @description Do a maximum likelihood Weibull mixture fit for the input data vector x using the Expectation-Maximization (EM) alogrithm.
#'
# @details
#'
#' @param x Locations at which to evaluate probability density function
#' @param z0 Initial value of mixture proportions
#' @param theta0 The value of shape and scale parameters with the ordering [sh1,sc1,sh2,sc2,...]
#' @param numSteps Number of EM steps to perform (default 25)
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return A list consisting of the final z (mixture proportions), theta (mixture parameters), and nll (negative log-likelihood).
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{calcPdfMatrixWeibMix,calcPdfWeibMix,calcNegLogLikWeibMix}
# @references

fitWeibMixMaxLikEM <- function(x, z0, theta0, numSteps = 25) {
  zLast <- z0 # Value of mixture proportions for preceding step
  thetaLast <- theta0 # Value of mixture parameters for preceding step
  nllLast <- calcNegLogLikWeibMix(theta0, x, z0) # Negative log-likelihood for preceding step

  for (ii in 1:numSteps) {
    # 1. Expectation update
    f <- calcPdfMatrixWeibMix(x, thetaLast)
    zNew <- colMeans(f / rowSums(f))

    # 2. Optimization step
    fit <- optim(thetaLast, calcNegLogLikWeibMix, x = x, z = zNew)
    thetaNew <- fit$par

    zLast <- zNew
    thetaLast <- thetaNew
    nllNew <- calcNegLogLikWeibMix(thetaNew, x, zNew)
    nllLast <- nllNew
  }
  return(list(z = zNew, theta = thetaNew, nll = nllNew))
}
