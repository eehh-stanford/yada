#' @title Calculate the posterior probability for the mixed cumulative probit model
#'
#' @description \code{calc_x_posterior} calculates the posterior probability (x|y,theta_x,theta_y) for the mixed cumulative probit model at the points in xcalc given the input parameterization y. xcalc is assumed to be evenly spaced.
#'
# @details
#'
#' @param xcalc A vector of ages at which to calculate the posterior probability
#' @param y An observation of traits used for the posterior calculation
#' @param theta_x Parameterization for prior on x
#' @param theta_y Parameterization for likelihood
#' @param hp Hyperparameters
#'
#'
#' @return A vector of posterior probabilities
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
calc_x_posterior <- function(xcalc,y,theta_x,theta_y,hp) {
  dx <- xcalc[2] - xcalc[1] # Assumes xcalc is evenly spaced
  p_xy <- calc_joint(xcalc,y,theta_x,theta_y,hp)
  p_x <- p_xy / sum(p_xy) / dx
  return(p_x)
}
