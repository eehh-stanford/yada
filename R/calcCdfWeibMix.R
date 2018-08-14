#' @title Calculate CDF given Weibull mixture proportions and parameters
#'  
#' @description Given the input vector x and Weibull mixture proportions and parameters z and theta, calculate the cumulative density function (CDF) for the input vector x.
#' 
# @details
#'
#' @param x Locations at which to evaluate cumulative density function
#' @param z Vector of mixture proportions
#' @param theta The value of shape and scale parameters with the ordering [sh1,sc1,sh2,sc2,...]
#' 
# @keywords
#' @export
#' 
# @examples
#' 
#' @return The CDF vector for the input vector x
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso \code{fitWeibMixMaxLikEM,calcPdfWeibMix,calcPdfMatrixWeibMix}
# @references 

calcCdfWeibMix <- function(x,z,theta) {
  numObs <- length(x) # Number of observations
  numMix <- length(theta)/2 # Number of mixtures
  # For x, the shape parameter, and scale parameter create a long vector of
  # length numObs*numMix for vectorized input to pweibull
  x_vect <- rep(x,numMix)
  shape_vect <- as.vector(t(matrix(rep(theta[seq(1,length(theta),by=2)],length(x)),nrow=numMix)))
  scale_vect <- as.vector(t(matrix(rep(theta[seq(2,length(theta),by=2)],length(x)),nrow=numMix)))
  # Calculate the cumulative density and transform to a matrix with dimensions numObs x numMix
  cdfMatrix <- matrix(pweibull(x_vect,shape_vect,scale_vect),ncol=numMix)
  return(rowSums(t(t(cdfMatrix)*z)))
}
