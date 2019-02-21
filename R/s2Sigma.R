#' @title Convert covariance matrix from vector to matrix representation
#'
#' @description \code{s2Sigma} converts the input vector s into an equivalent covariance matrix Sigma. The ordering of s is by row: Sigma_11, Sigma_12,...,Sigma_22,....
#'
#' @param s The covariance matrix in vector form
#' @param Sigma The covariance matrix in matrix form
#'
# @keywords
#' @export
#'
#' @return The covariance matrix Sigma
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

s2Sigma <- function(s) {
  
  K <- .5 * (sqrt(1+8*length(s))-1) # Number of variables
  Sigma <- matrix(rep(NA,K*K),nrow=K)
  Sigma[lower.tri(Sigma,diag=T)] <- s
  corr <- Sigma[lower.tri(Sigma)]
  Sigma <- t(Sigma)
  Sigma[lower.tri(Sigma)] <- corr
  Sigma <- t(Sigma)
  return(Sigma)
}
