#' @title Calculate the density for x given parameterization theta_x
#'
#' @description \code{calc_x_density} calculates the density at x given a parameterization for the density of theta_x. Currently, the exponential distribution, Weibull mixtures, and uniform distribution are supported
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
  } else if (tolower(theta_x$fitType) == 'offsetweibmix') { 
    return(calcPdfWeibMix(x+theta_x$weibOffset,theta_x$fit$lambda,c(rbind(theta_x$fit$shape,theta_x$fit$scale))))
  } else if (tolower(theta_x$fitType) == 'uniform') { 
    x <- rep(1,length(x))
    ind0 <- x < theta_x$xmin | theta_x$xmax < x
    x[ind0] <- 0
    return(x)
  } else {
    stop(paste('Unsupported fit type:',theta_x$fitType))
  }
}
