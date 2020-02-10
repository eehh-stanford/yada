#' @title Power law with a constant offset
#'
#' @description \code{powLaw} calculates the mean (h). \code{powLawSigma} calculates the noise (sigma, or sig for short). \code{powLawDensity} calculates the density. \code{powLawNegLogLik} calculates the negative log-likelihood. \code{fitPowLaw} returns the maximum likelihood fit. \code{simPowLaw} creates simulated data. 
#'
#' @details We assume that the response variable w is distributed as
#'
#' \deqn{w ~ N(h(x),sig(x)^2)}
#'
#' where x is the independent variable, h(x) the mean, sig(x) the noise, and N denotes a normal distribution. If sig is independent of x, the model is homoskedastic. Otherwise, it is heteroskedastic. For an observation (w,x), the likelihood is
#'
#' \deqn{l_w = 1/sqrt(2*pi)/sig*exp(-0.5*(w-h)^2/sig^2)}
#'
#' The negative log-likelihood is
#'
#' \deqn{eta_w = 0.5*log(2*pi) + log(sig) + 0.5*(w-h)^2/sig^2}
#'
#' For the mean and noise, we adopt the parametric forms
#'
#' \deqn{h = a*x^r + b}
#'
#' and
#'
#' \deqn{sig = s*(1 + kap*x)}
#' 
#' @param x Vector of independent variable observations
#' @param w Vector of dependent variable observations
#' @param a Multiplicative coefficient
#' @param r Scaling exponent
#' @param b Offset
#' @param s Baseline noise
#' @param kap Slope of noise (short for kappa) [Optional]
#' @param th_w Vector of parameters with ordering [a,r,b,s,kap]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLaw <- function(x,th_w) {
  # th_w has ordering [a,r,b,s,kap]
  return(th_w[1]*x^th_w[2] + th_w[3])
}

#' @export
powLawSigma <- function(x,th_w,hetero=F) {
  # th_w has ordering [a,r,b,s,kap]
  # returns a scalar for hetero=F even if x is not length 1
  sig <- th_w[4]
  if(hetero) {
    sig <- sig * (1+th_w[5]*x)
  }
  return(sig)
}

#' @export
powLawDensity <- function(x,w,th_w,hetero=F) {
  # th_w has ordering [a,r,b,s,kap]
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetero)
  1/sqrt(2*pi)/sig*exp(-0.5*(w-h)^2/sig^2)
}

#' @export
powLawNegLogLik <- function(th_w,x,w,hetero=F) {
  # th_w has ordering [a,r,b,s,kap]
  # eta_w is the negative log-likelihood
  # For optimization, th_w is the first input

  N <- length(x) # No error checking is done on input lengths
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetero)
  eta_w <- 0.5*log(2*pi)*N + sum(log(sig) + 0.5*(w-h)^2/sig^2)
  return(eta_w)
}


#' @export
fitPowLaw <- function(x,w,hetero=F) {
  # th_w has ordering [a,r,b,s,kap]

  # Initialize parameters
  r0 <- 1 # linear in x
  b0 <- min(w)
  a0 <- diff(range(w))/diff(range(x))
  s0 <- diff(range(w))/2

  th_w0 <- c(a0,r0,b0,s0)
  if(hetero) {
    gam0 <- 0
    th_w0 <- c(th_w0,gam0)
  }

  optimControl <- list(reltol=1e-12,maxit=100000)
  fit <- optim(th_w0,powLawNegLogLik,method='BFGS',control=optimControl,x=x,w=w,hetero=hetero,hessian=T)

  return(fit$par)
}

#' @export
simPowLaw <- function(N,th_x,th_w,hetero=F) {
  # N is the number of simulated observations
  # th_x parameterizes x. Currently, only a uniform distribution is supported
  # th_w parameterizes w (given x)
  a <- th_w[1]
  r <- th_w[2]
  b <- th_w[3]
  s <- th_w[4]

  x <- runif(N,th_x[1],th_x[2])
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetero)

  w <- h + rnorm(N)*sig
  
  return(list(x=x,w=w))
}


