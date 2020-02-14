#' @title Power law with ordinal observations
#'
#' @description \code{powLawOrd} calculates the mean (h). \code{powLawOrdSigma} calculates the noise (sigma, or sig for short). \code{powLawOrdNegLogLik} calculates the negative log-likelihood. \code{fitPowLawOrd} returns the maximum likelihood fit. \code{simPowLawOrd} creates simulated data. \code{powLawOrdCalc_x_list} transforms from a vector to list representation for the input data.
#'
#' @details We assume that the latent response variable vstar (for v^*) is distributed as
#'
#' \deqn{vstar ~ N(g(x),sig(x)^2)}
#'
#' where x is the independent variable, g(x) the mean, sig(x) the noise, and N denotes a normal distribution. If sig is independent of x, the model is homoskedastic. Otherwise, it is heteroskedastic. What is observed is not the latent response vstar, but rather an ordinal response variable v that is derived from vstar per
#'
#' \deqn{-Inf  < vstar <= tau_1 --> m = 0}
#' \deqn{tau_1 < vstar <= tau_2 --> m = 1}
#' \deqn{...}
#' \deqn{tau_M < vstar <= Inf   --> m = M}
#'
#' That is, there are M+1 ordinal categories m = 0,1,...M, and the ordered vector tau = [tau_1,...,tau_M] provides the cut-offs to derive the ordinal response, v, from the latent continuous response, vstar. Given this model, the likelihood of the observation (x,v) is
#'
#' \deqn{l_v = Phi((tau_{m+1} - g)/sig) - Phi((tau_m - g)/sig)}
#'
#' where Phi is the cumulative distribution function of the standard univariate normal with a mean of 0 and standard deviation of 1. We adopt the convention tau_0 = -Inf and tau_{M+1} = Inf and adopt the notation that tau_m is tau_lo and tau_{m+1} is tau_hi (and similarly for Phi and other variables). The negative log likelihood of the observation (x,v) is
#'
#' \deqn{eta_v = -log(Phi((tau_hi - g)/sig) - Phi((tau_lo - g)/sig))}
#' \deqn{      = -log(Phi_hi-Phi_lo)                                }
#'
#' For the mean and noise, we adopt the parametric forms
#'
#' \deqn{h = x^rho}
#'
#' and
#'
#' \deqn{sig = s*(1 + kap*x)}
#'
#' The choice of x^rho comes from using an offset power law, alpha*x^rho + beta, and requiring alpha=1 and beta=0 for identifiability. For optimization, it is often preferable to work with an unconstrained variable. This is supported via the optional input transformVar. rho, s, kap, and the differences between successive values of tau must be positive. This is accomplished by using, for example, rho_bar = log(rho) and rho = exp(rho_bar).
#'
#' @param x Vector of independent variable observations
#' @param vstar Vector of latent dependent variable observations (for v^*)
#' @param v Vector of dependent variable observations (ordinal)
#' @param rho Scaling exponent
#' @param s Baseline noise
#' @param kap Slope of noise (short for kappa) [Optional]
#' @param th_v Vector of parameters with ordering [rho,tau_1,...,tau_M,s,kap]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#' @param transformVar Whether a transformation of the parameterization is needed [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLawOrd <- function(x,th_v,transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  if(transformVar) {
    rho <- exp(th_v[1])
  } else {
    rho <- th_v[1]
  }
  return(x^rho)
}

#' @export
powLawOrdSigma <- function(x,th_v,hetero=F,transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # returns a scalar for hetero=F even if x is not length 1
  numParam <- length(th_v)

  if(hetero) {
    if(transformVar) { 
      s   <- exp(th_v[numParam-1])
      kap <- exp(th_v[numParam])
    } else {
      s   <- th_v[numParam-1]
      kap <- th_v[numParam]
    }
    sig <- s*(1 + kap*x)
  } else {
    if(transformVar) { 
      s   <- exp(th_v[numParam])
    } else {
      s   <- th_v[numParam]
    }
    sig <- s
  }

  return(sig)
}

#' @export
powLawOrdNegLogLik <- function(th_v,x_list,hetero=F,transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # eta_v is the negative log-likelihood
  # For optimization, make th_v the first input

  if(transformVar) {
    # Build hp
    if(!hetero) {
      hp <- list(paramModel='powLawOrdHomo')
      hp$J <- 1
      hp$M <- length(th_v) - 2
    } else {
      hp <- list(paramModel='powLawOrdHetero')
      hp$J <- 1
      hp$M <- length(th_v) - 3
    }
    th_v <- theta_y_unconstr2constr(th_v,hp)
  }

  M <- length(x_list) - 1 # number of ordinal categories is M+1
  rho <- th_v[1]         # rho
  tau <- th_v[2:(M+1)]   # tau_1 ... tau_M
  s   <- th_v[M+2]       # s

  if(hetero) {
    kap <- th_v[M+3]     # kappa
  }

  eta_v <- 0
 
  for(m in 0:M) { 
    x <- x_list[[m+1]]
    if(hetero) {
      sig <- s*(1+kap*x)
    } else {
      sig <- s 
    }

    if(m==0) {
      Phi_lo <- 0
    } else {
      Phi_lo <- pnorm( (tau[m]- x^rho)/sig )
    }

    if(m==M) {
      Phi_hi <- 1
    } else {
      Phi_hi <- pnorm( (tau[m+1]- x^rho)/sig )
    }

    eta_v <- eta_v - sum(log(Phi_hi - Phi_lo))
  }

  return(eta_v)
}

#' @export
fitPowLawOrd <- function(x,v,hetero=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  M <- length(unique(v)) - 1
  x_list <- powLawOrdCalc_x_list(x,v)

  if(hetero) {
    hp <- list(paramModel = 'powLawOrdHetero')
  } else {
    hp <- list(paramModel = 'powLawOrdHomo')
  }

  hp$J <- 1
  hp$M <- M
  # To initialize, fit a set of binary probits at the different levels in order
  # to get tau0 and s0

  # initial transition points (tau)
  tau0 <- rep(NA,M)
  s0Vect <- rep(NA,M)

  # Suppress the following warning (it's OK that it happens):
  # glm.fit: fitted probabilities numerically 0 or 1 occurred
  oldw <- getOption("warn")
  options(warn = -1)
  for(m in 0:(M-1)) {
    vbinary <- rep(0,length(v))
    vbinary[v > m] <- 1
    linProbit <- glm(vbinary~x,family=binomial(link='probit'))
    b <- as.numeric(linProbit$coefficients[1]) # intercept
    a <- as.numeric(linProbit$coefficients[2]) # slope
    tau0[m+1] <- -b/a
    s0Vect <-  1/a
  }
  options(warn = oldw)
  s0 <- mean(s0Vect)

  if(is.unsorted(tau0)) {
    warning('tau0 is not sorted. Attempting to resolve by sorting')
    tau0 <- sort(tau0)
  }

  rho0   <- 1

  th_v0 <- c(rho0,tau0,s0)
  if(hetero) {
    kap0   <- 0.001
    th_v0 <- c(th_v0,kap0)
  }

  th_v_bar0 <- theta_y_constr2unconstr(th_v0,hp)

  optimControl <- list(reltol=1e-12,maxit=100000)
  fit <- optim(th_v_bar0,powLawOrdNegLogLik,control=optimControl,x_list=x_list,hetero=hetero,hessian=T,transformVar=T,method='BFGS')
  th_v <- theta_y_unconstr2constr(fit$par,hp)
  
  return(list(fit=fit,th_v=th_v,th_v0=th_v0))
}

#' @export
simPowLawOrd <- function(N,th_x,th_v,hetero=F) {
  # N is the number of simulated observations
  # th_x parameterizes x. Currently, only a uniform distribution is supported
  # th_v parameterizes v (given x)

  if(hetero) {
    M <- length(th_v) - 3
  } else {
    M <- length(th_v) - 2
  }
  tau <- th_v[2:(M+1)]

  x <- runif(N,th_x[1],th_x[2])
  g   <- powLawOrd(x,th_v)
  sig <- powLawOrdSigma(x,th_v,hetero)

  vstar <- g + sig*rnorm(N)
  v <- rep(NA,N)

  for(n in 1:N) {
    v[n] <- as.numeric(cut(vstar[n],c(-Inf,tau,Inf))) - 1 # The ordinal observation
  }

  return(list(x=x,v=v,vstar=vstar))
}

#' @export
powLawOrdCalc_x_list <- function(x,v) {
  # For the input vectors x / v, create a list of length M + 1 in which the
  # entries of the list are the x-values for each ordinal category.
  M <- length(unique(v)) - 1 # Number of ordinal categories is M + 1
  x_list <- list()
  for(m in 0:M) {
    x_list[[m+1]] <- x[v == m]
  }

  return(x_list)
}
