#' @title Convert from covariance matrix to correlation vector
#'
#' @description \code{covMat2corrVect} converts the input covariance matrix Sigma into an equivalent vector of correlations (off diagonal terms). The diagonal terms should be $1$. The ordering is by row: Sigma_12, Sigma_13,...,Sigma_23,....
#'
#' @param Sigma The covariance matrix
#' @param s The correlation vector
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The correlation vector s
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

covMat2CorrVect <- function(Sigma) {
  return(Sigma[which(lower.tri(Sigma))])
}
