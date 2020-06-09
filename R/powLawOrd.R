#' @title Power law with ordinal observations
#'
#' @description \code{powLawOrd} calculates the mean (h). \code{powLawOrdSigma} calculates the noise (sigma, or sig for short). \code{powLawOrdNegLogLikVect} calculates a vector of the negative log-likelihoods. \code{powLawOrdNegLogLik} calculates the negative log-likelihood (sum of \code{powLawOrdNegLogLikVect}). \code{fitPowLawOrd} returns the maximum likelihood fit. \code{simPowLawOrd} creates simulated data. \code{powLawOrdCalc_x_list} transforms from a vector to list representation for the input data.
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
#' \deqn{sig = s*(1 + kappa*x)}
#'
#' The choice of x^rho comes from using an offset power law, alpha*x^rho + beta, and requiring alpha=1 and beta=0 for identifiability. For optimization, it is often preferable to work with an unconstrained variable. This is supported via the optional input transformVar. rho, s, kappa, and the differences between successive values of tau must be positive. This is accomplished by using, for example, rho_bar = log(rho) and rho = exp(rho_bar).
#'
#' @param x Vector of independent variable observations
#' @param vstar Vector of latent dependent variable observations (for v^*)
#' @param v Vector of dependent variable observations (ordinal)
#' @param rho Scaling exponent
#' @param s Baseline noise
#' @param kappa Slope of noise [Optional]
#' @param th_v Vector of parameters with ordering [rho,tau_1,...,tau_M,s,kappa]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#' @param transformVar Whether a transformation of the parameterization is needed [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLawOrd <- function(x,th_v,transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kappa]
  if(transformVar) {
    rho <- exp(th_v[1])
  } else {
    rho <- th_v[1]
  }
  return(x^rho)
}

#' @export
powLawOrdSigma <- function(x,th_v,hetSpec='none',transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kappa,lambda]
  # returns a scalar for hetero=F even if x is not length 1
  hetero <- hetSpec != 'none'
  numParam <- length(th_v)

  if(hetero) {
    if(hetSpec == 'sd_pow') {
      if(transformVar) { 
        s   <- exp(th_v[numParam-2])
        kappa <- exp(th_v[numParam-1])
        lambda <- exp(th_v[numParam])
      } else {
        s   <- th_v[numParam-2]
        kappa <- th_v[numParam-1]
        lambda <- th_v[numParam]
      }
    } else {
      if(transformVar) { 
        s   <- exp(th_v[numParam-1])
        kappa <- exp(th_v[numParam])
      } else {
        s   <- th_v[numParam-1]
        kappa <- th_v[numParam]
      }
    }

    if(hetSpec == 'sd_x') {
      sig <- s*(1+kappa*x)
    } else if(hetSpec == 'sd_resp') {
      sig <- s*(1+kappa*x^rho)
    } else if(hetSpec == 'sd_pow') {
      sig <- s*(1+kappa*x^lambda)
    } else {
      stop(paste('Unrecognized hetSpec,',hetSpec))
    }
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
powLawOrdNegLogLikVect <- function(th_v,x,v,hetSpec='none',transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kappa]
  # eta_v is the negative log-likelihood
  # For optimization, make th_v the first input

  hetero <- hetSpec != 'none'

  M <- calc_M(th_v,hetSpec)

  if(transformVar) {
    # Build modSpec
    modSpec <- list(meanSpec='powLaw')
    modSpec$J <- 1
    modSpec$M <- M
    modSpec$hetSpec <- hetSpec
    if(hetero) {
      modSpec$hetGroups <- 1
    }

    th_v <- theta_y_unconstr2constr(th_v,modSpec)
  }

  rho <- th_v[1]           # rho
  tau <- th_v[2:(M+1)]     # tau_1 ... tau_M
  s   <- th_v[M+2]         # s

  if(hetero) {
    kappa <- th_v[M+3]     # kappa
    if(hetSpec == 'sd_pow') {
      lambda <- th_v[M+4]  # lambda
    }
  }

  # Initialize vector of negative log-likelihoods (eta_v)
  N <- length(x)
  if(N != length(v)) {
    stop('Input vectors x and v do not match in length')
  }
  eta_v <- rep(NA,N)
 
  for(m in 0:M) { 
    ind <- v == m
    if(sum(ind) > 0 ) {
      xm <- x[ind]
      if(hetero) {
        if(hetSpec == 'sd_x') {
          sig <- s*(1+kappa*xm)
        } else if(hetSpec == 'sd_resp') {
          sig <- s*(1+kappa*xm^rho)
        } else if(hetSpec == 'sd_pow') {
          sig <- s*(1+kappa*xm^lambda)
        } else {
          stop(paste('Unrecognized hetSpec,',hetSpec))
        }
      } else {
        sig <- s 
      }

      if(m==0) {
        Phi_lo <- 0
      } else {
        Phi_lo <- pnorm( (tau[m]- xm^rho)/sig )
      }

      if(m==M) {
        Phi_hi <- 1
      } else {
        Phi_hi <- pnorm( (tau[m+1]- xm^rho)/sig )
      }

      eta_v[ind] <- - log(Phi_hi - Phi_lo)
    }
  }

  return(eta_v)
}

#' @export
powLawOrdNegLogLik <- function(th_v,x,v,hetSpec='none',transformVar=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kappa]
  # eta_v is the negative log-likelihood
  # For optimization, make th_v the first input
  return(sum(powLawOrdNegLogLikVect(th_v,x,v,hetSpec,transformVar)))
}

#' @export
calc_M <- function(th_v,hetSpec) {
  # A helper function to calculate M from the length of th_v
  hetero <- hetSpec != 'none'
  if(hetero) {
    if(hetSpec == 'sd_pow') {
      M <- length(th_v) - 4
    } else { 
      M <- length(th_v) - 3
    }
  } else {
    M <- length(th_v) - 2
  }
  return(M)
}

#' @export
fitPowLawOrd <- function(x,v,hetSpec='none') {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  hetero <- hetSpec != 'none'
  M <- length(unique(v)) - 1

  # Build modSpec
  modSpec <- list(meanSpec='powLaw')
  modSpec$J <- 1
  modSpec$M <- M
  modSpec$hetSpec <- hetSpec
  if(hetero) {
    modSpec$hetGroups <- 1
  }

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
    kappa0   <- 0.001
    th_v0 <- c(th_v0,kappa0)
    if(hetSpec == 'sd_pow') {
      lambda0   <- 1
      th_v0 <- c(th_v0,lambda0)
    }
  }

  th_v_bar0 <- theta_y_constr2unconstr(th_v0,modSpec)

  optimControl <- list(reltol=1e-12,maxit=100000,ndeps=rep(1e-8,length(th_v_bar0)))
  fit <- optim(th_v_bar0,powLawOrdNegLogLik,control=optimControl,x=x,v=v,hetSpec=hetSpec,hessian=T,transformVar=T,method='BFGS')

  th_v <- theta_y_unconstr2constr(fit$par,modSpec)
  
  return(th_v)
}

#' @export
simPowLawOrd <- function(N,th_x,th_v,hetSpec='none') {
  # N is the number of simulated observations
  # th_x parameterizes x. Currently, only a uniform distribution is supported
  # th_v parameterizes v (given x)

  hetero <- hetSpec != 'none'

  M <- calc_M(th_v,hetSpec)
  tau <- th_v[2:(M+1)]

  x <- runif(N,th_x[1],th_x[2])
  g   <- powLawOrd(x,th_v)
  sig <- powLawOrdSigma(x,th_v,hetSpec)

  vstar <- g + sig*rnorm(N)
  v <- rep(NA,N)

  for(n in 1:N) {
    v[n] <- as.numeric(cut(vstar[n],c(-Inf,tau,Inf))) - 1 # The ordinal observation
  }

  return(list(x=x,v=v,vstar=vstar))
}
