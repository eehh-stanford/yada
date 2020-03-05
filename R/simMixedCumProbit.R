#' @title Created simulated data for a mixed cumulative probit analysis
#'
#' @description \code{simMixedCumProbit} creates simulated data for a mixed cumulative probit analysis.
#'
#' @param theta_y_list theta_y as a list
#' @param N number of simulated observations
#' @param theta_x Rate parameter for independent variable draw from exponential distribution
#' @param hp Hyperparameters (needed to guide generation of simulated data)
#' @param missDataProp Proportion of missing data [default 0]
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
  if(theta_x$fitType == 'uniform') {
    x <- runif(N,theta_x$xmin,theta_x$xmax)
  } else {
    stop('Only uniform currently supported')
  }

  hetero <- is_hetero(hp$paramModel)

  # Create latent observations
  Ystar <- matrix(NA,hp$J+hp$K,N)
  for(n in 1:N) {
    covMat <- theta_y_list$Sigma
    if(hetero) {
      covMat <- covMat * (1 + theta_y_list$kappa*x[n])
    }
    Ystar[,n] <- t(MASS::mvrnorm(1,mu=rep(0,hp$J+hp$K),Sigma=covMat))
  }
  # add offset
  for(rr in 1:nrow(Ystar)) {
    if(rr <= hp$J) {
      j <- rr
      rho <- theta_y_list$rho[j]
      Ystar[rr,] <- Ystar[rr,] + x[n]^rho
    } else {
      k <- rr - hp$J
      a <- theta_y_list$a[k]
      r <- theta_y_list$r[k]
      b <- theta_y_list$b[k]
      Ystar[rr,] <- Ystar[rr,] + a*x[n]^r + b
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
