#' @title Maximum likelihood fit for Siler distribution
#'
#' @description Do a maximum likelihood Siler fit for the input data vector x. To allow unconstrained optimization to be used, optimization is done on the transformed variable abar = log(a/a0).
#'
# @details
#'
#' @param x Locations at which to evaluate probability density function
#' @param a0 Initial value for optimization. Default from Gage and Dyke 1986, Table 2, Level 15
#' @param x0 x-offset. The default is 0 (no offset).
#' @param calcHessian Whether to calculate the Hessian
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return A list consisting of the fit (on the transformed variable abar) and maximum likelihood estimate of a. Optionally, the Hessian of the log-likelihood is returned, which allows estimation of the standard errors of the maximum likelihood estimate via the observed Fisher information matrix.
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{dsiler}
# @references

fitSilerMaxLik <- function(x, a0 = c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1), x0 = 0, calcHessian = FALSE,hessEps=1e-10) {
  nllbar <- function(abar, x, a0) {
    a <- a0 * exp(abar)
    -sum(log(dsiler(x, a, x0)))
  }

  fit <- optim(rep(0, 5), nllbar, x = x, a0 = a0)

  afit <- a0 * exp(fit$par)
  if (calcHessian) {
    ll <- function(a) { # log likelihood as a function of a
      sum(log(dsiler(x, a, x0)))
    }
    H <- numDeriv::hessian(ll, afit,method.args=list(eps=hessEps))
    return(list(fit = fit, a = afit, hessian = H))
  } else {
    return(list(fit = fit, a = afit))
  }
}
