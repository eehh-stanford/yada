#' @title Calculate posterior for mixed cumulative probit model
#'
#' @description \code{calc_mixed_cum_probit_posterior} calculates the posterior probability of an observation as a function of a vector xv, which is assumed to be uniformly sampled.
#'
# @details
#'
#' @param theta_y theta_y as a list (single theta_y sample) or list of lists (multiple theta_y samples)
#' @param theta_x Parameterization for distribution of x (exponential rate)
#' @param xv A vector of ages at which to calculate posterior
#' @param y A new observation
#' @param hp Hyperparameters
#' @param asLog A logical variable indicating whether to return the (unnormalized) log of the likelihood
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
calc_mixed_cum_probit_posterior <- function(theta_y,theta_x,xv,y,hp,asLog=F) {
  '%dopar%' <- foreach::'%dopar%'

  # Determine if there is one or multiple theta_y
  singleVal <- 'Sigma' %in% names(theta_y)

  if(singleVal) {
    logPost <- foreach::foreach(n=1:length(xv), .combine=cbind) %dopar% {
      logPost <- calcLogLik_theta_y(theta_y,xv[n],y,hp)
      logPost <- logPost + dexp(xv[n],theta_x,log=T)
    }
    logPost <- as.numeric(logPost)
    if(!asLog) {
      dx <- xv[2] - xv[1]
      post <- exp(logPost)
      post <- post/dx/sum(post)
    } else {
      post <- logPost
    }
  } else {
    logPostMat <- matrix(NA,length(theta_y),length(xv))
    for(tt in 1:length(theta_y)) {
      logPostMat[tt,] <- calc_mixed_cum_probit_posterior(theta_y[[tt]],theta_x,xv,y,hp,asLog=T)
    }
    post <- rep(NA,length(xv))
    for(n in 1:length(xv)) {
      post[n] <- matrixStats::logSumExp(logPostMat[,n]) - log(length(theta_y))
    }
    if(!asLog) {
      post <- exp(post)
      dx <- xv[2] - xv[1]
      post <- post/dx/sum(post)
    }
  }
  return(post)
}
