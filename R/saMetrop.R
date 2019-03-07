#' @title Metropolis sampling for adaptive simulated annealing algorithm
#'
#' @description This function samples using the Metropolis algorithm using a fixed input covariance matrix for the proposal distribution. The distribution being sampled from is set by the input cost function and temperature.
#'
# @details
#'
#' The first input to the target distribution function, targetDist, must be the variable to be sampled. Additional inputs to targetDist can be input as named variables.
#'
#' The `control` argument is a list that can supply any of the following (otherwise defaults, in brackets, are be used):
#' \itemize{
#'   \item{`numSamp`}{Number of samples after initialization [1000]}
#'   \item{`verbose`}{Whether to print out information as the run proceeds [F]}
#'   \item{`fileName`}{Filename for saving}
#'   \item{`savePeriod`}{Period (of samples) for saving [1000]}
#' }
#'
#' @param costFunc The cost function (often a negative log-likelihood)
#' @param X_0 The starting point for sampling
#' @param ... Further arguments to be passed to costFunc
#' @param control A list of control parameters. See Details.
#'
#' @return A list containing the samples along with summary information
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
saMetrop <- function(costFunc,X_0,temp,...,control=list()) {

  # If necessary, set control defaults
  if(!('numSamp' %in% names(control))) {
    control$numSamp <- 1000
  }

  if(!('verbose' %in% names(control))) {
    control$verbose <- F
  }

  if(!('C' %in% names(control))) {
    if(length(X_0) == 1) {
      control$C <- 1e-6
    } else {
      control$C <- diag(c(rep(1e-6,length(X_0))))
    }
  }

  if(!('temp' %in% names(control))) {
    control$temp <- 1
  }

  cf <- function(X) costFunc(X,...)

  # Initialize variables
  X_t <- X_0
  cost_t <- cf(X_t)

  X_mat <- matrix(NA,length(X_0),control$numSamp)
  costVect <- vector()
  acceptVect <- vector()
  for(tt in 1:control$numSamp) {
    if(control$verbose) {
      print('----')
      print(tt)
      print(temp)
    }

    if(length(X_t)==1) {
      X_tp1 <- X_t +  rnorm(1,sd=sqrt(control$C))
    } else {
      X_tp1 <- MASS::mvrnorm(1,mu=X_t,Sigma=control$C)
    }

    
    # Calculate the acceptance ratio, alpha
    # By construction, the proposal distribution is symmetric. Hence:
    # A = p(theta_tp1|alpha,D) / p(theta_t|alpha,D)
    # A = p(D|theta_tp1) * p(theta_tp1|alpha) / 
    #     p(D|theta_t) / p(theta_t|alpha) 
    cost_tp1 <- cf(X_tp1)
    if(!is.finite(cost_tp1)) {
      accept <- F
    } else {
      alpha <- min(1,exp(-(cost_tp1-cost_t)/temp))
      accept <- runif(1) < alpha
    }

    acceptVect[tt] <- accept
    if(!accept) {
      X_tp1 <- X_t
      cost_tp1 <- cost_t
    }
    X_mat[,tt] <- X_tp1
    costVect[tt] <- cost_tp1

    # Get ready for next sample
    X_t <- X_tp1
    cost_t <- cost_tp1
    if(control$verbose) {
      print(cost_t)
      print(accept)
    }
  }
  return(list(X_mat=X_mat,cost=costVect,acceptVect=acceptVect))
}
