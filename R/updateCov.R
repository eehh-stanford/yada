#' @title Recursively update a covariance matrix
#'
#' @description Recursively update a covariance matrix (and mean) by adding new observations.
#'
#' @param X A new measurement
#' @param covObj A list with the current covariance matrix, mean, and number of observations
#'
#' @return covObj, a list with the updated covariance matrix, mean, and number of observations
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
updateCov <- function(X,covObj=NA) {

  if(all(is.na(covObj))) {
    covObj <- list(mean=X,cov=NA,n=1)
    return(covObj)
  }

  covObj$n <- covObj$n + 1 # Update number of observations

  if(covObj$n==2) {
    X1 <- covObj$mean
    covObj$mean <- X1/2 + X/2
    dX1 <- X1 - covObj$mean
    dX2 <- X - covObj$mean
    covObj$cov <- tcrossprod(dX1,dX1) + tcrossprod(dX2,dX2)
    return(covObj)
  }

  dx <- covObj$mean - X # previous mean minus new X
  covObj$cov <- covObj$cov * (covObj$n-2)/(covObj$n-1) + tcrossprod(dx,dx)/covObj$n
  covObj$mean <- covObj$mean*(covObj$n-1)/covObj$n + X/covObj$n
  return(covObj)
}
