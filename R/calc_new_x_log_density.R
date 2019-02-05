#' @title Calculate the log density for a new observation x
#'
#' @description \code{calc_new_x_log_density} calculates the log density for a new observation x given known y. This is in essence a wrapper for MCMC sampling of x for the probability density p(x|y,D,alpha). Optionally, the input can be the transformed variable xbar = log(x). The densities of x and x_bar are realated by p_xbar(xbar) = p_x(xbar)*exp(xbar).
#'
# @details
#'
#' @param xbar log of the scalar age at which to calculate the log density
#' @param y The known value of the response variables
#' @param theta_y theta_y as a list (single theta_y sample) or list of lists (multiple theta_y samples)
#' @param theta_x Parameterization for distribution of x (exponential rate)
#' @param hp Hyperparameters
#' @param transformVar Boolean indicating whether the input is xbar [optional; default FALSE]
#'
# @keywords
#' @export
#'
#' @return The log density
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
calc_new_x_log_density <- function(x,y,theta_y,theta_x,hp,transformVar=F) {
  if(!transformVar) {
    if(x < 0) {
      return(-Inf)
    } else {
      return(calc_mixed_cum_probit_posterior(theta_y,theta_x,x,y,hp,asLog=T))
    }
  } else {
    return(calc_mixed_cum_probit_posterior(theta_y,theta_x,exp(x),y,hp,asLog=T) + x)
  }
}
