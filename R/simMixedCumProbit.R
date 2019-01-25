#' @title Created simulated data for a mixed cumulative probit analysis
#'
#' @description \code{simMixedCumProbit} creates simulated data for a mixed cumulative probit analysis.
#'
#' @param theta_y_list theta_y as a list
#' @param N number of simulated observations
#' @param theta_x Rate parameter for independent variable draw from exponential distribution
#' @param hp Hyperparameters (needed to guide generation of simulated data)
#' @param missDataProp Proportion of missing data [default 0]

#' @param hp
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return sim A list of simulated data variables
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

simMixedCumProbit <- function(theta_y_list,N,theta_x,hp,missDataProp=0) {
  # Draw for the independent variable
  # scalar x variable
  x <- rexp(N,theta_x)

  # Create latent observations
  Ystar <- t(MASS::mvrnorm(N,mu=rep(0,hp$J+hp$K),Sigma=theta_y_list$Sigma))
  # add offset
  for(rr in 1:dim(Ystar)[1]) {
    if(rr <= hp$J) {
      param <- c(theta_y_list$alpha[rr],theta_y_list$rho[rr])
      Ystar[rr,] <- Ystar[rr,] + genCrraOrd(x,param)
    } else {
      param <- c(theta_y_list$a[rr-hp$J],theta_y_list$r[rr-hp$J],theta_y_list$b[rr-hp$J])
      Ystar[rr,] <- Ystar[rr,] + genCrra(x,param)
    }
  }

  # actual observations without missing data
  Y0 <- Ystar

  for(j in 1:hp$J) {
    for(n in 1:N) {
      Y0[j,n] <- as.numeric(cut(Ystar[j,n],c(-Inf,theta_y_list$tau[[j]],Inf))) - 1 # The ordinal observation
    }
  }

  if(missDataProp > 0) {
    # Add missing data
    Y <- Y0
    Y[runif(nrow(Y0)*ncol(Y0),0,1) <= missDataProp] <- NA
    return(list(x=x,Y=Y,Ystar=Ystar,Y0=Y0))
  } else {
    return(list(x=x,Ystar=Ystar,Y=Y0))
  }

}
