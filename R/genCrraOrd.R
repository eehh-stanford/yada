#' @title Ordinal generalized CRRA parametric model
#'
#' @description \code{genCrraOrd} calculates the response (y). \code{genCrraOrdProb} calculates the category probabilities. fitGenCrraOrd returns the least squares fit. simGenCrraOrd creates simulated data.
#'
#'
#' @details The parametric form of the generalized Coefficient of Relative Risk Aversion (CRRA) for ordinal data is (no slope or intercept):
#'
#' \deqn{y = x^(1-\rho)}
#'
#' where \deqn{x} is the independent variable, \deqn{y} the latent dependent (response) variable, \deqn{\alpha} the multiplicative coefficient, \deqn{\rho} the scaling exponent, and  there is no offset.
#'
#' @param x Vector of independent variable observations
#' @param y Vector of dependent variable observations
#' @param \deqn{\rho} Scaling exponent
#' @param \deqn{\sigma} Standard deviation of noise term [sig in code]
#' @param \deqn{\tau} Vector of category boundaries (one fewer elements than number of categories)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @export
genCrraOrd <- function(x,param,transformVar=F) {
  # param has ordering [rho,sigma,tau]
  if(!transformVar) {
    rho <- param[1]
  } else {
    rho <- 1 - exp(-param[1])
  }
  return(x^(1-rho))
}

#' @export
genCrraOrdProb <- function(param,x,y,transformVar=F) {
  # param has ordering [rho,sig,tau]
  J <- length(param) - 2
  sig <- param[2]
  if(!transformVar) {
    rho <- param[1]
    tau <- param[3:length(param)]
  } else {
    rho <- 1 - exp(-param[1])
    tau0 <- param[3]
    dtau <- exp(param[4:length(param)])
    tau <- tau0 + c(0,cumsum(dtau))
  }
  M <- length(tau)
  lo <- rep(-Inf,length(x))
  lo[y != 0] <- tau[y[y!=0]]
  hi <- rep(Inf,length(x))
  hi[y != J] <- tau[y[y!=J]+1]
  return(pnorm(hi,genCrraOrd(x,rho,transformVar=F),sig) - pnorm(lo,genCrraOrd(x,rho,transformVar=F),sig))
}

#' @export
genCrraOrdLik <- function(param,x,y,transformVar=F) {
  # param has ordering [rho,sig,tau]
  return(sum(log(genCrraOrdProb(param,x,y,transformVar))))
}

#' @export
fitGenCrraOrd <- function(x,y) {
  # Remove possible missing value
  B <- is.na(x) | is.na(y)
  x <- x[!B]
  y <- y[!B]

  M <- max(y) # Categories are m = 0,1,...,M
  paramCont <- fitGenCrraNoIntercept(x,y)
  xm <- rep(NA,M+1)
  ym <- rep(NA,M+1)
  yPred <- list()
  tauVect <- rep(NA,M)
  for(m in 0:M) {
    xm[m+1] <- mean(x[y == m])
    ym[m+1] <- xm[m+1]^paramCont[2]
    yPred[[m+1]] <- genCrra(x[y==m],c(1,paramCont[2],0))
    if(m > 0) {
      tauVect[m] <- mean(c(mean(yPred[[m]]),mean(yPred[[m+1]])))
    }
  }
  ycat <- rep(NA,length(x))
  for(n in 1:length(x)) {
    ycat[n] <- as.numeric(cut(genCrra(x[n],c(1,paramCont[2],0)),c(-Inf,tauVect,Inf))) - 1 # The ordinal observation
  }

  propRight <- sum(ycat == y) / length(y)
  tauScale <- mean(diff(tauVect))/2 # Characteristic tau scale for random draws
  param0 <- c(-log(1-paramCont[2]),-tauScale / qnorm(propRight/2,0,1),tauVect[1],log(diff(tauVect)))

  #renorm <- max(y-genCrra(x,paramCont))
  #param0 <- c(paramCont[1]/renorm,paramCont[2],1,rep(log(1/renorm),M-1))
  #param0 <- c(1/mean(x),0,1,rep(0,M-1))

  fit <- optim(param0,genCrraOrdLik,x=x,y=y,transformVar=T,control=list(fnscale=-1))
  param <- fit$par
  param[1] <- 1 - exp(-param[1])
  tau0 <- param[3]
  dtau <- exp(param[4:length(param)])
  tau <- tau0 + c(0,cumsum(dtau))
  param[3:length(param)] <- tau

  #ycat2 <- rep(NA,length(x))
  #for(n in 1:length(x)) {
  #  ycat2[n] <- as.numeric(cut(genCrra(x[n],c(1,param[2],0)),c(-Inf,param[3:length(param)],Inf))) - 1 # The ordinal observation
  #}
  #propRight2 <- sum(ycat2 == y) / length(y)
  return(param)
}

#' @export
simGenCrraOrd <- function(x,param,latent=F) {
  ystar <- genCrraOrd(x,param) + param[2]*rnorm(length(x))
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


#' @export
genCrraNoInterceptRmse <- function(param,x,y,transformW=F) {
  # param has ordering [a,r]
  return(sqrt(mean((y - genCrra(x,c(param[1],param[2],0),transformW))^2)))
}

#' @export
fitGenCrraNoIntercept <- function(x,y) {
  # Remove possible missing value
  B <- is.na(x) | is.na(y)
  x <- x[!B]
  y <- y[!B]

  # Calculate a starting parameterization to get rough scale
  indMin <- which.min(x)
  indMax <- which.max(x)
  a0 <- (y[indMax] - y[indMin]) / (x[indMax] - x[indMin])
  param0 <- c(a0,0)
  fit <- optim(param0,genCrraNoInterceptRmse,x=x,y=y,transformW=T)
  param <- as.numeric(c(fit$par[1],1-exp(-fit$par[2])))
  return(param)
}

##' @export
#genCrraOrdJoint <- function(x,y,fit_y,theta_x) {
#  p_x <- calc_x_density(x,theta_x)
#  p_ygx <- genCrraOrdProb(fit_y,x,y)
#  p_xy <- p_x * p_ygx
#  return(p_xy)
#}
