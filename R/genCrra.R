#' @title Generalized CRRA parametric model
#'
#' @description \code{genCrra} calculates the response (y). fitGenCrra returns the least squares fit. \code{simGenCrra} creates simulated data. 
#'
#'
#' @details The parametric form of the generalized Coefficient of Relative Risk Aversion (CRRA) is:
#'
#' \deqn{y = a*x^(1-r) + b}
#'
#' where \deqn{x} is the independent variable, \deqn{y} the dependent (response) variable, \deqn{a} the multiplicative coefficient, \deqn{r} the scaling exponent, and \deqn{b} the offset.
#'
#' @param x Vector of independent variable observations
#' @param y Vector of dependent variable observations
#' @param a Multiplicative coefficient
#' @param r Scaling exponent
#' @param b Offset
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @export
genCrra <- function(x,param,transformW=F) {
  # param has ordering [a,r,b]
  if(!transformW) {
    r <- param[2]
  } else {
    r <- 1 - exp(-param[2])
  }
  return(param[1] * x^(1-r) + param[3])
}

#' @export
genCrraRmse <- function(param,x,y,transformW=F,useSig=F) {
  # param has ordering [a,r,b]
  return(sqrt(mean((y - genCrra(x,param,transformW))^2)))
}

#' @export
genCrraLogLik <- function(param,x,y,transformW=F) {
  # param has ordering [a,r,b,sig]
  logLik <- -length(x)*(log(2*pi)/2 + log(param[4])) - sum((y - genCrra(x,param,transformW))^2) / 2 / param[4]^2
  return(logLik)
}

genCrraLik <- function(param,x,y,transformW=F,asLog=F) {
  # param has ordering [a,r,b,sig]
  lik <- dnorm(y,genCrra(x,param[1:3]),param[4],log=asLog)
  return(lik)
}

#' @export
fitGenCrra <- function(x,y,param0=NA,fitSig=F) {
  # Remove possible missing value
  B <- is.na(x) | is.na(y)
  x <- x[!B]
  y <- y[!B]

  if(all(is.na(param0))) {
    # Calculate a starting parameterization to get rough scale
    indMin <- which.min(x)
    indMax <- which.max(x)
    a0 <- (y[indMax] - y[indMin]) / (x[indMax] - x[indMin])
    b0 <- 0.5 * (y[indMin] + y[indMax] - a0 * (x[indMin] + x[indMax]) )
    param0 <- c(a0,0,b0,1)
  } else {
    param0 <- c(param0[1],-log(1-param0[2]),param0[3],param0[4])
    if(fitSig) {
      param0 <- c(param0,1)
    }
  }

  if(fitSig) {
    fit <- optim(param0[1:3],genCrraRmse,x=x,y=y,transformW=T)
  } else {
    fit <- optim(param0,genCrraRmse,x=x,y=y,transformW=T)
  }
  param <- as.numeric(c(fit$par[1],1-exp(-fit$par[2]),fit$par[3]))
  if(fitSig) {
    param <- c(param,sqrt(mean((y-genCrra(x,param))^2)))
  }
  
  return(param)
}

#' @export
simGenCrra <- function(x,param) {
  useSig <- length(param) > 3
  if(!useSig) {
    y <- genCrra(x,param) + rnorm(length(x))
  } else {
    y <- genCrra(x,param[1:3]) + param[4]*rnorm(length(x))
  }
  return(y)
}
