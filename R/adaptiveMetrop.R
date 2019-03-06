#' @title Adaptive Metropolis sampling
#'
#' @description This function implements the adaptive Metropolis algorithm described by Harrio, Saksman, and Tamminen (2001), An adaptive Metropolis algorithm, in Bernoulli, Volume 7, Number 2, pages 223 through 242. A starting point for the variable being sampled must be given, but otherwise the algorithm can sample with default settings.
#'
# @details
#'
#' The first input to the target distribution function, targetDist, must be the variable to be sampled. Additional inputs to targetDist can be input as named variables.
#'
#' The `control` argument is a list that can supply any of the following (otherwise defaults, in brackets, are be used):
#' \itemize{
#'   \item{`numInit`}{Number of samples using initial proposal distribution [100]}
#'   \item{`numSamp`}{Number of samples after initialization [1000]}
#'   \item{`verbose`}{Whether to print out information as the run proceeds [F]}
#'   \item{`fileName`}{Filename for saving}
#'   \item{`savePeriod`}{Period (of samples) for saving [1000]}
#' }
#'
#' @param targetDist The target distribution [pi in Harrio et al. (2001)]
#' @param X_0 The starting point for sampling
#' @param ... Further arguments to be passed to targetDist
#' @param control A list of control parameters. See Details.
#'
#' @return A list containing the samples along with summary information
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
adaptiveMetrop <- function(targetDist,X_0,...,control=list()) {

  # If necessary, set control defaults
  if(!('numInit' %in% names(control))) {
    control$numInit <- 100
  }

  if(!('numSamp' %in% names(control))) {
    control$numSamp <- 1000
  }

  if(!('verbose' %in% names(control))) {
    control$verbose <- F
  }

  if(!('C_0' %in% names(control))) {
    if(length(X_0) == 1) {
      control$C_0 <- 1e-6
    } else {
      control$C_0 <- diag(c(rep(1e-6,length(X_0))))
    }
  }

  if(!('s_d' %in% names(control))) {
    control$s_d <- (2.4)^2 / length(X_0)
  }

  if(!('epsilon' %in% names(control))) {
    control$epsilon <- 1e-12
  }

  if(!('temp' %in% names(control))) {
    control$temp <- 1
  }

  #if(!('scale' %in% names(control))) {
  #  control$scale <- 1
  #}

  I_d <- diag(length(X_0))
  td <- function(X) targetDist(X,...)

  # Initialize variables
  X_t <- X_0
  targetVal_t <- td(X_0)

  targetVal_0 <- targetVal_t

  meanX <- X_t

  X_list <- list()
  targetValVect <- vector()
  for(tt in 1:(control$numInit + control$numSamp)) {
    if(tt <= control$numInit) {
      if(control$verbose) {
        print('---- [burn in]')
      }
      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(control$C_0))
      } else {
        X_tp1 <- MASS::mvrnorm(1,mu=X_t,Sigma=control$C_0)
      }

    } else {
      if(control$verbose) {
        print('----')
      }
      C_t <- control$s_d * (covX + control$epsilon*I_d)
      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(C_t))
      } else {
        X_tp1 <- MASS::mvrnorm(1,mu=X_t,Sigma=C_t)
      }
    }
    if(control$verbose) {
      print(tt)
    }
    
    # Calculate the acceptance ratio, alpha
    # By construction, the proposal distribution is symmetric. Hence:
    # A = p(theta_tp1|alpha,D) / p(theta_t|alpha,D)
    # A = p(D|theta_tp1) * p(theta_tp1|alpha) / 
    #     p(D|theta_t) / p(theta_t|alpha) 
    targetVal_tp1 <- td(X_tp1)
    if(!is.finite(targetVal_tp1)) {
      accept <- F
    } else {
      alpha <- min(1,targetVal_tp1/targetVal_t)
      if(is.na(alpha)) {
        print(targetVal_t)
        print(targetVal_tp1)
      }
      accept <- runif(1) < alpha
    }

    if(!accept) {
      X_tp1 <- X_t
      targetVal_tp1 <- targetVal_t
    }
    X_list[[tt]] <- X_tp1
    targetValVect[tt] <- targetVal_tp1

    # Get ready for next sample
    X_t <- X_tp1
    targetVal_t <- targetVal_tp1
    meanX <- meanX * tt/(tt+1) + X_t/(tt+1)
    if(tt == 1) {
      B_t <- tcrossprod(X_0,X_0) + tcrossprod(X_t,X_t)
    } else {
      B_t <- B_t*(tt-1)/tt + tcrossprod(X_t,X_t)/tt
    }
    covX <- B_t - (tt+1)*tcrossprod(meanX,meanX)/tt
    #covX <- (covX + t(covX))/2 # ensures covX remains positive definite
    if(control$verbose) {
      print(log(targetVal_t)*control$temp)
      print(accept)
    }
  }
  return(list(X_list=X_list,target=targetValVect,C_t=C_t))
}
