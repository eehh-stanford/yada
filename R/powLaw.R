#' @title Power law with a constant offset
#'
#' @description \code{powLaw} calculates the mean (h). \code{powLawSigma} calculates the noise (sigma, or sig for short). \code{powLawDensity} calculates the density. \code{powLawNegLogLikVect} calculates a vector of negative log-likelihood. \code{powLawNegLogLik} calculates the negative log-likelihood (sum of \code{powLawNegLogLikVect}). \code{fitPowLaw} returns the maximum likelihood fit. \code{simPowLaw} creates simulated data. 
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
#' \deqn{sig = s*(1 + kappa*x)}
#' 
#' @param x Vector of independent variable observations
#' @param w Vector of dependent variable observations
#' @param a Multiplicative coefficient
#' @param r Scaling exponent
#' @param b Offset
#' @param s Baseline noise
#' @param kappa Slope of noise [Optional]
#' @param th_w Vector of parameters with ordering [a,r,b,s,kappa]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#' @param transformVar Whether a transformation of the parameterization is needed [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLaw <- function(x,th_w,transformVar=F) {
  # th_w has ordering [a,r,b,s,kappa]
  a <- th_w[1]
  r <- th_w[2]
  b <- th_w[3]
  if(transformVar) {
    a <- exp(a)
    r <- exp(b)
  }
  return(a*x^r + b)
}

#' @export
powLawSigma <- function(x,th_w,hetSpec='none',transformVar=F) {
  # th_w has ordering [a,r,b,s,kappa]
  # returns a scalar for hetero=F even if x is not length 1
  sig <- th_w[4]
  if(transformVar) {
    sig <- exp(sig)
  }

  if(!is_th_w_hetero(th_w)) {
    return(sig)
  }

  # If this point is reached, the model is heteroskedastic
  a <- th_w[1]
  r <- th_w[2]
  b <- th_w[3] 
  kappa <- th_w[5]

  if(transformVar) {
    a <- exp(a)
    r <- exp(r)
    kappa <- exp(kappa)
  }

  if(hetSpec == 'sd_x') {
    sig <- sig * (1+kappa*x)
  } else if(hetSpec == 'sd_resp') {
    sig <- sig * (1+a*kappa*x^r)
    sig <- sig * (1+a*kappa*x^l)
  } else {
    stop(paste('Unrecognized hetSpec,',hetSpec))
  }
  return(sig)
}

#' @export
powLawDensity <- function(x,w,th_w,hetSpec='none') {
  # th_w has ordering [a,r,b,s,kappa]
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetSpec)
  return(1/sqrt(2*pi)/sig*exp(-0.5*(w-h)^2/sig^2))
}

#' @export
powLawNegLogLikVect <- function(th_w,x,w,hetSpec='none',transformVar=F) {
  # th_w has ordering [a,r,b,s,kappa]
  # eta_w is the negative log-likelihood
  # For optimization, th_w is the first input
  hetero <- is_th_w_hetero(th_w)

  if(transformVar) {
    # Build modSpec
    modSpec <- list(meanSpec='powLaw')
    modSpec$K <- 1
    modSpec$hetSpec <- hetSpec
    if(hetero) {
      modSpec$hetGroups <- 1
    }

    th_w <- theta_y_unconstr2constr(th_w,modSpec)
  }

  N <- length(x) # No error checking is done on input lengths
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetSpec)
  eta_w <- 0.5*log(2*pi) + log(sig) + 0.5*(w-h)^2/sig^2
  return(eta_w)
}

#' @export
powLawNegLogLik <- function(th_w,x,w,hetSpec='none',transformVar=F) {
  # th_w has ordering [a,r,b,s,kappa]
  # eta_w is the negative log-likelihood
  # For optimization, th_w is the first input
  return(sum(powLawNegLogLikVect(th_w,x,w,hetSpec,transformVar)))
}

#' @export
fitPowLaw <- function(x,w,hetSpec='none') {
  # th_w has ordering [a,r,b,s,kappa]
  hetero <- hetSpec != 'none'

  # Initialize parameters
  r0 <- 1 # linear in x
  b0 <- min(w)
  a0 <- diff(range(w))/diff(range(x))
  s0 <- diff(range(w))/2

  th_w0 <- c(a0,r0,b0,s0)
  if(hetero) {
    kappa0 <- 0.0001
    th_w0 <- c(th_w0,kappa0)
  }

 # Build modSpec
  modSpec <- list(meanSpec='powLaw')
  modSpec$K <- 1
  modSpec$hetSpec <- hetSpec
  if(hetero) {
    modSpec$hetGroups <- 1
  }

  th_w_bar0 <- theta_y_constr2unconstr(th_w0,modSpec)
  optimControl <- list(reltol=1e-12,maxit=100000,ndeps=rep(1e-8,length(th_w_bar0)))
  fit <- optim(th_w_bar0,powLawNegLogLik,method='BFGS',control=optimControl,x=x,w=w,hessian=T,hetSpec=hetSpec,transformVar=T)

  th_w <- theta_y_unconstr2constr(fit$par,modSpec)

  return(th_w)
}

#' @export
simPowLaw <- function(N,th_x,th_w,hetSpec='none') {
  # N is the number of simulated observations
  # th_x parameterizes x. Currently, only a uniform distribution is supported
  # th_w parameterizes w (given x)
  a <- th_w[1]
  r <- th_w[2]
  b <- th_w[3]
  s <- th_w[4]

  hetero <- hetSpec != 'none'
  if(hetero) {
    kappa <- th_w[5]
  }

  x <- runif(N,th_x[1],th_x[2])
  h   <- powLaw(x,th_w)
  sig <- powLawSigma(x,th_w,hetSpec)

  w <- h + rnorm(N)*sig
  
  return(list(x=x,w=w))
}

#' @export
is_th_w_hetero <- function(th_w) {
  if(length(th_w) == 4) {
    return(F)
  } else if(length(th_w) == 5) {
    return(T)
  } else {
    stop(paste('Length of th_w should be 4 or 5, not',length(th_w)))
  }
}
