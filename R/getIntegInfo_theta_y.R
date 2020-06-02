#' @title Get integration info for log-likelihood calculation for theta_y
#'
#' @description \code{getIntegInfo_theta_y} determines the integration limits for the log-likelihood calculation for theta_y
#'
#' @details
#' Determine the integration limits. For ordinal variables, this comes from tau in theta. For missing data, the integration limit is -Inf to Inf. Two variables are returned in a list:
#'
#' (a) doIntegral, which is a [(J+K) x N] matrix J+K is the number of variables and N is the number of observations. doIntegral is T if an integration must be done and F otherwise. F occurs for non-missing, continuous observations.
#'
#' (b) intLim, which is a [(J+K x N x 2] array and provides the integration limits for the maximum likelihood calculation
#'
#' @param theta_y_list theta_y as a list
#' @param Y Matrix of dependent variable observations
#'
# @keywords
#' @export
#'
#' @return A list consisting of the matrix doIntegral and array intLim [see details]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

getIntegInfo_theta_y <- function(theta_y_list,Y) {

  J <- get_J(theta_y_list$modSpec)
  K <- get_K(theta_y_list$modSpec)
  haveOrd  <- J != 0
  haveCont <- K != 0
  numVar <- J+K

  isCont <- c(rep(F,J),rep(T,K)) # Is the variable continuous?
  if(is.matrix(Y)) {
    # If a matrix is input, it should have dimensions (J+K) x N
    if(nrow(Y) != (J+K)) {
      stop('Number of rows in Y should equal J+K')
    }
    N <- ncol(Y) # number of observations
  } else {
    # The input is a vector (or a scalar)
    if(length(Y) != (J+K)) {
      stop('Input vector (or scalar) length should equal J+K')
    }
    N <- 1
  }
  # For non-missing observations, columns of doIntegral are T for ordinal variables and F for continuous variables
  # doIntegral has dimensions (J+K) x N
  if(N > 1) {
    doIntegral <- matrix(rep(!isCont,N),nrow=numVar)
  } else {
    doIntegral <- !isCont
  }
  # The other possibility for a T is that there is missing data
  doIntegral[is.na(Y)] <- T
  
  # Create the limit array for observations requiring integration
  # limArray has dimensions (J+K) x N x 2
  if(N > 1) {
    limArray <- array(rep(NA,numVar*N*2),dim=c(numVar,N,2))
    # Add limits for missing variables
    limArray[,,1][is.na(Y)] <- -Inf
    limArray[,,2][is.na(Y)] <-  Inf
  } else {
    limArray <- matrix(NA,length(Y),2)
    # Add limits for missing variables
    limArray[,1][is.na(Y)] <- -Inf
    limArray[,2][is.na(Y)] <-  Inf
  }


  # Add limits for non-missing ordinal variables
  if(haveOrd) {
    for(j in 1:J) {
      if(J > 1) {
        tau <- c(-Inf,theta_y_list$tau[[j]],Inf)
      } else {
        tau <- c(-Inf,unlist(theta_y_list$tau),Inf)
      }
      if(N > 1) {
        for(n in 1:N) {
          if(!is.na(Y[j,n])) {
            limArray[j,n,1] <- tau[Y[j,n]+1]
            limArray[j,n,2] <- tau[Y[j,n]+2]
          }
        }
      } else {
          if(!is.na(Y[j])) {
            limArray[j,1] <- tau[Y[j]+1]
            limArray[j,2] <- tau[Y[j]+2]
          }
      }
    }
  }
  return(list(doIntegral=doIntegral,limArray=limArray))
}
