#' @title Plot posterior for mixed cumulative probit model
#'
#' @description \code{plot_mixed_cum_probit_posterior} calculates and plots the posterior probability of an observation as a function of a vector xv, which is assumed to be uniformly sampled.
#'
# @details
#'
#' @param theta_y_list theta_y as a list
#' @param theta_x Parameterization for distribution of x (exponential rate)
#' @param xv A vector of ages at which to calculate posterior
#' @param y A new observation
#' @param hp Hyperparameters
#' @param knownX Known x value [optional]
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
plot_mixed_cum_probit_posterior <- function(theta_y_list,theta_x,xv,y,hp,knownX=NA) {
  fv <- calc_mixed_cum_probit_posterior(theta_y_list,theta_x,xv,y,hp)
  xMin <- min(xv)
  xMax <- max(xv)
  yMin <- 0
  yMax <- max(fv)
  if(!is.na(knownX)) {
    xMax <- max(xMax,knownX)
  }

  plot(1,type="n",xlab="x", ylab='density',xlim=c(xMin,xMax),ylim=c(yMin,yMax))
  if(!is.na(knownX)) {
    lines(c(1,1)*knownX,c(0,yMax),col='grey',lwd=3)
  }
  lines(xv,fv,lwd=3)
}
