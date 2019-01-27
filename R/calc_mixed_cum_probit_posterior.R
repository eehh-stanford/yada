#' @title Calculate posterior for mixed cumulative probit model
#'
#' @description \code{calc_mixed_cum_probit_posterior} calculates the posterior probability of an observation as a function of a vector xv, which is assumed to be uniformly sampled.
#'
# @details
#'
#' @param theta_y_list theta_y as a list
#' @param theta_x Parameterization for distribution of x (exponential rate)
#' @param xv A vector of ages at which to calculate posterior
#' @param y A new observation
#' @param hp Hyperparameters
#'
# @keywords
#' @export
#'
#' @return The vector of posterior probabilities
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
# @references

#' @export
calc_mixed_cum_probit_posterior <- function(theta_y_list,theta_x,xv,y,hp) {
   logPost <- foreach::foreach(n=1:length(xv), .combine=cbind) %dopar% {
    logPost <- calcLogLik_theta_y(theta_y_list,xv[n],y,hp)
    logPost <- logPost + dexp(xv[n],theta_x)
  }
  dx <- xv[2] - xv[1]
  post <- exp(logPost)
  post <- post/dx/sum(post)
  return(post)
}
