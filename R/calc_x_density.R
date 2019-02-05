#' @title Calculate the density for x given parameterization theta_x
#'
#' @description \code{calc_x_density} calculates the density at x given a parameterization for the density of theta_x. Currently, the exponential distribution and Weibull mixtures are supported
#'
# @details
#'
#' @param x A vector of ages at which to calculate the density
#' @param theta_x List with the fit type and parameter vector
#'
#'
#' @return The vector of densities
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
# @references

#' @export
calc_x_density <- function(x,theta_x) {
  if(tolower(theta_x$fitType) == 'exponential') { 
    return(dexp(x,theta_x))
  } else if (tolower(theta_x$fitType) == 'weibmix') { 
    return(calcPdfWeibMix(x,theta_x$fit$z,theta_x$fit$theta))
  } else {
    stop(paste('Unsupported fit type:',theta_x$fitType))
  }
}
