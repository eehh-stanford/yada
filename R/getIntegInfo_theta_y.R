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
  J <- length(theta_y_list$alpha)
  K <- length(theta_y_list$a)
  isCont <- c(rep(F,J),rep(T,J)) # Is the variable continuous?
  # For non-missing observations, columns of doIntegral are T for ordinal variables and F for continuous variables
  # doIntegral has dimensions (J+K) x N
  if(is.matrix(Y)) {
    doIntegral <- matrix(rep(!isCont,ncol(Y)),nrow=nrow(Y))
  } else {
    doIntegral <- !isCont
  }
  # The other possibility for a T is that there is missing data
  doIntegral[is.na(Y)] <- T
  
  # Create the limit array for observations requiring integration
  # limArray has dimensions (J+K) x N x 2
  if(is.matrix(Y)) {
    limArray <- array(rep(NA,nrow(Y)*ncol(Y)*2),dim=c(nrow(Y),ncol(Y),2))
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
  for(j in 1:J) {
    tau <- c(-Inf,theta_y_list$tau[[j]],Inf)
    if(is.matrix(Y)) {
      for(n in 1:ncol(Y)) {
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
  return(list(doIntegral=doIntegral,limArray=limArray))
}
