#' @title Convert from correlation vector to covariance matrix
#'
#' @description \code{corrVect2CovMat} converts the input correlation vector into an equivalent covariance matrix Sigma. The diagonal terms of Sigma are $1$. The ordering is by row: Sigma_12, Sigma_13,...,Sigma_23,....
#'
#' @param s The correlation vector
#' @param Sigma The covariance matrix
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The covariance matrix Sigma
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

corrVect2CovMat <- function(s) {
  K <- .5 * (1 + sqrt(1+8*length(s))) # Number of variables
  Sigma <- diag(rep(1,K))
  counter <- 1
  # iterate over rows less than K
  for(rr in 1:(K-1)) {
      # iterate over columns greater than the row number, rr
      for(cc in (rr+1):K) {
          Sigma[rr,cc] <- s[counter]
          Sigma[cc,rr] <- s[counter]
          counter <- counter + 1
      }
  }
  return(Sigma)
}
