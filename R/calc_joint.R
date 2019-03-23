#' @title Calculate the joint probability for the mixed cumulative probit model
#'
#' @description \code{calc_joint} calculates the joint density probability p(x,y|theta_x,theta_y) for the mixed cumulative probit model at the points in xcalc. xcalc is assumed to be evenly spaced.
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
#' @return A vector of joint probabilities
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
calc_joint <- function(xcalc,y,theta_x,theta_y,hp) {
  logPrior <- log(calc_x_density(xcalc,theta_x))

  # Since y is a single observation, use a parallel for loop
  '%dopar%' <- foreach::'%dopar%'
  logLik <- foreach::foreach(n=1:length(xcalc), .combine=cbind) %dopar% {
    ll <- calcLogLik_theta_y(theta_y,xcalc[n],y,hp)
  }
  logJoint <- logLik + logPrior 
  #for(n in 1:length(xcalc)) {
  #  logJoint[n] <- logJoint[n] + calcLogLik_theta_y(theta_y,xcalc[n],y,hp)
  #}

  fv <- exp(logJoint)
  fv[!is.finite(fv)] <- 0
  return(fv)
}
