#' @title Convert theta_y from list to vector format
#'
#' @description \code{theta_yList2Vect} converts the variable theta_y formatted as a list to a vector. The ordering of the vector is [rho, a, r, b, tau, s].
#'
#' @param theta_y_list theta_y as a list
#' @param theta_y theta_y as a vector
#'
# @keywords
#' @export
#'
#' @return theta_y as a vector
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

theta_yList2Vect <- function(theta_y_list) {
  haveOrd <- 'rho' %in% names(theta_y_list)
  haveCont <- 'r' %in% names(theta_y_list)

  if(haveOrd) {
    tau <- unlist(theta_y_list$tau)
  }
  s <- Sigma2s(theta_y_list$Sigma)

  if(haveOrd && haveCont) {
    return(c(theta_y_list$rho,theta_y_list$a,theta_y_list$r,theta_y_list$b,tau,s,theta_y_list$gamma))
  } else if(haveOrd && !haveCont) {
    return(c(theta_y_list$rho,tau,s))
  } else {
    # !haveOrd && haveCont
    return(c(theta_y_list$a,theta_y_list$r,theta_y_list$b,s,theta_y_list$gamma))
  }
}
