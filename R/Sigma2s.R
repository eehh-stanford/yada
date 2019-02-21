#' @title Convert covariance matrix from matrix to vector representation
#'
#' @description \code{Sigma2s} converts the input covariance matrix Sigma into an equivalent vector s. The ordering of s is by row: Sigma_11, Sigma_12,...,Sigma_22,....
#'
#' @param Sigma The covariance matrix in matrix form
#' @param s The covariance matrix in vector form
#'
# @keywords
#' @export
#'
#' @return The vector s representing the covariance matrix  
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

Sigma2s <- function(Sigma) {
  
  return(Sigma[lower.tri(Sigma,diag=T)])
}
