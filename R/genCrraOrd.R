#' @title Ordinal generalized CRRA parametric model
#'
#' @description \code{genCrraOrd} calculates the response (y). \code{genCrraOrdProb} calculates the category probabilities. fitGenCrraOrd returns the least squares fit. simGenCrraOrd creates simulated data.
#'
#'
#' @details The parametric form of the generalized Coefficient of Relative Risk Aversion (CRRA) for ordinal data is (no intercept):
#'
#' \deqn{y = \alpha*x^(1-\rho)}
#'
#' where \deqn{x} is the independent variable, \deqn{y} the latent dependent (response) variable, \deqn{\alpha} the multiplicative coefficient, \deqn{\rho} the scaling exponent, and  there is no offset.
#'
#' @param x Vector of independent variable observations
#' @param y Vector of dependent variable observations
#' @param \deqn{\alpha} Multiplicative coefficient
#' @param \deqn{\rho} Scaling exponent
#' @param \deqn{\tau} Vector of category boundaries (one fewer elements than number of categories)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @export
genCrraOrd <- function(x,param,transformVar=F) {
  # param has ordering [alpha,rho,tau]
  if(!transformVar) {
    rho <- param[2]
  } else {
    rho <- 1 - exp(-param[2])
  }
  return(param[1] * x^(1-rho))
}

#' @export
genCrraOrdProb <- function(param,x,y,transformVar=F) {
  # param has ordering [alpha,rho,tau]
  J <- length(param) - 2
  alpha <- param[1]
  if(!transformVar) {
    rho <- param[2]
    tau <- param[3:length(param)]
  } else {
    rho <- 1 - exp(-param[2])
    tau0 <- param[3]
    dtau <- exp(param[4:length(param)])
    tau <- tau0 + c(0,cumsum(dtau))
  }
  M <- length(tau)
  lo <- rep(-Inf,length(x))
  lo[y != 0] <- tau[y[y!=0]]
  hi <- rep(Inf,length(x))
  hi[y != J] <- tau[y[y!=J]+1]
  return(pnorm(hi,genCrraOrd(x,c(alpha,rho),transformVar=F)) - pnorm(lo,genCrraOrd(x,c(alpha,rho),transformVar=F)))
}

#' @export
genCrraOrdLik <- function(param,x,y,transformVar=F) {
  # param has ordering [alpha,rho,tau]
  return(sum(log(genCrraOrdProb(param,x,y,transformVar))))
}

#' @export
fitGenCrraOrd <- function(x,y) {
  # Remove possible missing value
  B <- is.na(x) | is.na(y)
  x <- x[!B]
  y <- y[!B]

  M <- max(y) # Categories are m = 0,1,...,M
  param0 <- c(1/mean(x),0,1,rep(0,M-1))

  fit <- optim(param0,genCrraOrdLik,x=x,y=y,transformVar=T,control=list(fnscale=-1))
  param <- fit$par
  param[2] <- 1 - exp(-param[2])
  tau0 <- param[3]
  dtau <- exp(param[4:length(param)])
  tau <- tau0 + c(0,cumsum(dtau))
  param[3:length(param)] <- tau
  return(param)
}

#' @export
simGenCrraOrd <- function(x,param,latent=F) {
  ystar <- genCrraOrd(x,param) + rnorm(length(x))
  N <- length(x)
  y <- rep(NA,N)
  tau <- param[3:length(param)]
  for(n in 1:N) {
    y[n] <- as.numeric(cut(ystar[n],c(-Inf,tau,Inf))) - 1 # The ordinal observation
  }
  if(latent) {
    return(list(y=y,ystar=ystar))
  } else {
    return(y)
  }
}
