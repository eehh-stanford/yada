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
saMetrop <- function(costFunc,init,temp=NA,...,control=list()) {

  # init is either a vector (X_0) or list (continuing a chain)
  haveChain <- is.list(init)
  if(!haveChain) {
    X_0 <- init
  } else {
    X_0 <- as.vector(init$X_mat[,ncol(init$X_mat)])
    control0 <- control
    control <- init$control # These will get reset if the input control is not empty
  }

  # If necessary, set control defaults
  if(is.na(temp)) {
    if(!haveChain) {
      stop('Cannot initialize a new chain without an input temperature')
    }
    temp <- init$temp
  }

  if(!haveChain && !('numSamp' %in% names(control))) {
    control$numSamp <- 1000
  }

  if(haveChain && ('numSamp' %in% names(control0))) {
    control$numSamp <- control0$numSamp
  }

  if(!haveChain && !('verbose' %in% names(control))) {
    control$verbose <- F
  }

  if(haveChain && ('verbose' %in% names(control0))) {
    control$verbose <- control0$verbose
  }


  if(!haveChain && !('C_0' %in% names(control))) {
    if(length(X_0) == 1) {
      control$C_0 <- 1e-6
    } else {
      control$C_0 <- diag(c(rep(1e-6,length(X_0))))
    }
  }

  if(haveChain && ('C_0' %in% names(control0))) {
    control$C_0 <- control0$C_0
  }

  if(!haveChain && !('t0' %in% names(control))) {
    control$t0 <- 100
  }

  if(haveChain && ('t0' %in% names(control0))) {
    control$t0 <- control0$t0
  }

  if(!haveChain && !('s_d' %in% names(control))) {
    control$s_d <- (2.4)^2 / length(X_0)
  }

  if(haveChain && ('s_d' %in% names(control0))) {
    control$s_d <- control0$s_d
  }

  if(!haveChain && !('epsilon' %in% names(control))) {
    control$epsilon <- 1e-12
  }

  if(haveChain && ('epsilon' %in% names(control0))) {
    control$epsilon <- control0$epsilon
  }

  I_d <- diag(length(X_0))
  cf <- function(X) costFunc(X,...)

  # Initialize variables
  X_t <- X_0
  cost_t <- cf(X_t)

  if(!haveChain) {
    covObj <- updateCov(X_t)
    ttOffset <- 0
  } else {
    covObj <- init$covObj
    ttOffset <- ncol(init$X_mat)
  }

  X_mat <- matrix(NA,length(X_0),control$numSamp)
  costVect <- vector()
  acceptVect <- vector()
  for(tt in 1:control$numSamp) {
    if(tt+ttOffset <= control$t0) {    
      if(control$verbose) {
        print('-- C_0 --')
      }

      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(control$C_0))
      } else {
        X_tp1 <- MASS::mvrnorm(1,mu=X_t,Sigma=control$C_0)
      }
    } else {
      # tt > t0 [accounting for offset]
      if(control$verbose) {
        print('-- C_t --')
      }

      C_t <- control$s_d * (covObj$cov + control$epsilon*I_d)
      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(C_t))
      } else {
        X_tp1 <- MASS::mvrnorm(1,mu=X_t,Sigma=C_t)
      }
    }
    if(control$verbose) {
      print(tt+ttOffset)
      print(temp)
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
    covObj <- updateCov(X_t,covObj)
    if(control$verbose) {
      print(cost_t)
      print(accept)
    }
  } # end main loop
  if(haveChain) {
    X_mat <- cbind(init$X_mat,X_mat)
    costVect <- c(init$costVect,costVect)
    acceptVect <- c(init$acceptVect,acceptVect)
  }
  returnList <- list(X_mat=X_mat,costVect=costVect,acceptVect=acceptVect,covObj=covObj,control=control,temp=temp)
  if(!haveChain) {
    returnList$firstX <- X_0
  } else {
    returnList$firstX <- init$firstX
  }
  return(returnList)
}
