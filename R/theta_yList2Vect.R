#' @title Convert theta_y from list to vector format
#'
#' @description \code{theta_yList2Vect} converts the variable theta_y formatted as a list to a vector. The ordering of the vector is [alpha, rho, a, r, b, tau, s].
#'
#' @param theta_y_list theta_y as a list
#' @param theta_y theta_y as a vector
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return theta_y as a vector
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

theta_yList2Vect <- function(theta_y_list) {
  tau <- unlist(theta_y_list$tau)
  s <- covMat2CorrVect(theta_y_list$Sigma)
  return(c(theta_y_list$alpha,theta_y_list$rho,theta_y_list$a,theta_y_list$r,theta_y_list$b,tau,s))
}
