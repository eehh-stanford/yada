#' @title Sample theta_y in the mixed cumulative probit model
#'
#' @description \code{sample_theta_y} uses the adaptive Metropolis algorithm of Haario et al. 2001, https://projecteuclid.org/euclid.bj/1080222083, to sample the vector theta_y for the cumulative probit model.
#'
# @details
#'
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param hp Hyperparameters
#' @param numSamp Number of samples
#' @param burnIn Number of burn-in samples
#' @param verbose Whether to print out information as the run proceeds [optional]
#' @param start Starting value of theta_y_list for sampling [optional]
#' @param fileName File name for saving [optional]
#' @param savePeriod Period (of samples) for saving [optiona; default 1000]
#' @param known Known values, a list with the field theta_y, which is also a list [optional]
#' @param varNames Variable names for saving to file [optional]
#'
#' @return theta_y as a list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
estimate_theta_y <- function(x,Y,hp,verbose=T) {
  # Initialize theta_y and get ready for sampling
  theta_y0_list <- init_theta_y(x,Y,hp$J)
  #theta_y0 <- theta_yList2Vect(theta_y0_list)
  phi_y0_list <- theta_y_list2phi_y_list(theta_y0_list)
  phi_y0 <- phi_y_list2phi_y(phi_y0_list)
  
  # set up bounds
  lower <- c()
  upper <- c()
  # rho
  lower <- c(lower,rep(-10,hp$J))
  upper <- c(upper,rep(1,hp$J))

  # a
  lower <- c(lower,-10*rep(max(abs(phi_y0_list$a)),hp$K))
  upper <- c(upper,10*rep(max(abs(phi_y0_list$a)),hp$K))

  # r
  lower <- c(lower,-10*rep(max(abs(phi_y0_list$r)),hp$K))
  upper <- c(upper,10*rep(max(abs(phi_y0_list$r)),hp$K))

  # b
  lower <- c(lower,-10*rep(max(abs(phi_y0_list$b)),hp$K))
  upper <- c(upper,10*rep(max(abs(phi_y0_list$b)),hp$K))

  # tau_bar
  tau0Max <- -Inf
  dtauMax <- -Inf
  for(j in 1:hp$J) {
    tau_bar <- phi_y0_list$tau_bar[[j]]
    tau0Max <- max(tau0Max,abs(tau_bar[1]))
    dtauMax <- max(dtauMax,abs(tau_bar[2:length(tau_bar)]))
  }
  low2 <- c()
  upp2 <- c()
  for(j in 1:hp$J) {
    low2 <- c(low2,-10*tau0Max,-10*rep(dtauMax,hp$M[j]-1))
    upp2 <- c(upp2,10*tau0Max,10*rep(dtauMax,hp$M[j]-1))
  }
  lower <- c(lower,low2)
  upper <- c(upper,upp2)

  Umax <- max(abs(c(phi_y0_list$psi,phi_y0_list$omega)))
  lower <- c(lower,rep(0,length(phi_y0_list$psi)))
  upper <- c(upper,rep(10*Umax,length(phi_y0_list$psi)))
  lower <- c(lower,rep(-10*Umax,length(phi_y0_list$omega)))
  upper <- c(upper,rep(10*Umax,length(phi_y0_list$omega)))

  fit <- GenSA(lower=lower,upper=upper,fn=phi_yNegLogLikWrapper,x=x,Y=Y,hp=hp)
  return(theta_y)
}

#' @export
phi_yNegLogLikWrapper <- function(phi_y,x,Y,hp,checkValid=T) {
  theta_y_list <- phi_y2theta_y_list(phi_y,hp)
  if(!theta_yIsValid(theta_y_list,hp)) {
    return(Inf)
  }

  logLik <- calcLogLik_theta_y(theta_y_list,x,Y,hp)
  if(~is.finite(logLik)) {
    return(Inf)
  }
  return(-logLik)
}

#' @export
phi_y2theta_y_list <- function(phi_y,hp) {
 theta_y_list <- list(paramModel=hp$paramModel)
  
  last <- 0 # The last index added
  number <- hp$J # Number of elements to add

  # Add ordinal variable
  theta_y_list$rho <- phi_y[(last+1):(last+number)]

  # Add continuous variables
  last <- last + number
  number <- hp$K
  theta_y_list$a <- phi_y[(last+1):(last+number)]

  last <- last + number
  theta_y_list$r <- phi_y[(last+1):(last+number)]

  last <- last + number
  theta_y_list$b <- phi_y[(last+1):(last+number)]

  # Add tau
  tau <- list()
  for(j in 1:hp$J) {
    last <- last + number
    number <- hp$M[j]
    tau0 <- phi_y[(last+1)]
    dtau <- exp(phi_y[(last+2):(last+number)])
    tau[[j]] <- tau0 + c(0,cumsum(dtau))
  }
  theta_y_list$tau <- tau

  # Add covariance matrix
  last <- last + number
  number <- hp$K+hp$J
  psi <- phi_y[(last+1):(last+number)]

  last <- last + number
  number <- choose(hp$K+hp$J,2)
  omega <- phi_y[(last+1):(last+number)]
  U <- diag(psi)
  U[upper.tri(U)] <- omega
  theta_y_list$Sigma <- t(U) %*% U
  return(theta_y_list)
}


#' @export
theta_y_list2phi_y_list <- function(theta_y_list) {
  phi_y_list <- list(paramModel=theta_y_list$paramModel) # For regularized generalized CRRA
  phi_y_list$rho <- theta_y_list$rho
  phi_y_list$a <- theta_y_list$a
  phi_y_list$r <- theta_y_list$r
  phi_y_list$b <- theta_y_list$b

  phi_y_list$tau_bar <- list()
  for(j in 1:length(theta_y_list$tau)) {
    tau <- theta_y_list$tau[[j]]
    tau0 <- tau[1]
    dtau <- diff(tau)
    phi_y_list$tau_bar[[j]] <- c(tau0,log(dtau))
  }
  
  U <- chol(theta_y_list$Sigma)
  phi_y_list$psi <- diag(U)
  phi_y_list$omega <- U[upper.tri(U)]
  return(phi_y_list)
}

phi_y_list2phi_y <- function(phi_y_list) {
  tau_bar <- unlist(phi_y_list$tau_bar)
  return(c(phi_y_list$rho,phi_y_list$a,phi_y_list$r,phi_y_list$b,tau_bar,phi_y_list$psi,phi_y_list$omega))
}
