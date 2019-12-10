#' @title Power law with ordinal observations
#'
#' @description \code{powLaw} calculates the mean (h). \code{powLawSigma} calculates the noise (sigma, or sig for short). \code{powLawNegLogLik} calculates the negative log-likelihood. \code{powLawGradNegLogLik} calculates the gradient of the negative log-likelihood. \code{fitPowLaw} returns the maximum likelihood fit. \code{simPowLawOrd} creates simulated data. 
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
#' The choise of x^rho comes from using an offset power law, alpha*x^rho + beta, and requiring alpha=1 and beta=0 for identifiability. The components of the gradient of eta_v are
#'
#' \deqn{eta_v_rho    =  delta_phi / delta_Phi * x^rho * log(x)          }
#' \deqn{eta_v_tau_lo =  phi_lo / delta_Phi / sig                        }
#' \deqn{eta_v_tau_hi = -phi_hi / delta_Phi / sig                        }
#' \deqn{eta_v_s      =  1/delta_Phi *                                   }
#' \deqn{                   (phi_hi*(tau_hi - x^rho)                     }
#' \deqn{                  - phi_lo*(tau_lo - x^rho)) / sig^2 * (1+kap*x)}
#' \deqn{eta_v_kap    =  1/delta_Phi *                                   }
#' \deqn{                   (phi_hi*(tau_hi - x^rho)                     }
#' \deqn{                  - phi_lo*(tau_lo - x^rho)) / sig^2 * s * x    }
#'
#' where phi is the density function of the standard univariate normal, delta_Phi = Phi_hi - Phi_lo, and delta_phi = phi_hi - phi_lo.
#'
#' @param x Vector of independent variable observations
#' @param vstar Vector of latent dependent variable observations (for v^*)
#' @param v Vector of dependent variable observations (ordinal)
#' @param rho Scaling exponent
#' @param s Baseline noise
#' @param kap Slope of noise (short for kappa) [Optional]
#' @param th_v Vector of parameters with ordering [rho,tau_1,...,tau_M,s,kap]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLawOrd <- function(x,th_v) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  return(x^th_v[1])
}

#' @export
powLawOrdSigma <- function(x,th_v,hetero=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # returns a scalar for hetero=F even if x is not length 1
  numParam <- length(th_v)
  if(hetero) {
    sig <- th_v[numParam-1]*(1 + th_v[numParam]*x)
  } else {
    sig <- th_v[numParam]
  }

  return(sig)
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


#' @export
powLawOrdNegLogLik <- function(th_v,x_list,hetero=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # eta_v is the negative log-likelihood
  # For optimization, th_v is the first input

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
powLawOrdGradNegLogLik <- function(th_v,x_list,hetero=F) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # eta_v is the negative log-likelihood
  # For optimization, th_v is the first input

  M <- length(x_list) - 1 # number of ordinal categories is M+1
  rho <- th_v[1]         # rho
  tau <- th_v[2:(M+1)]   # tau_1 ... tau_M
  s   <- th_v[M+2]       # s
  if(hetero) {
    kap <- th_v[M+3]     # kappa
  }

 
  if(hetero) {
    gradVect <- rep(0,3+M)
  } else {
    gradVect <- rep(0,2+M)
  }
  
  for(m in 0:M) { 
    x <- x_list[[m+1]]
    x_to_rho <- x^rho

    # Handle the special case x = 0 by setting log_x to zero for these cases
    log_x <- log(x)
    log_x[x==0] <- 0

    if(hetero) {
      sig <- s*(1+kap*x)
    } else {
      sig <- s
    }

    sig_sq <- sig^2
   
    if(m == 0) {
      tau_lo <- 0 # only used when phi_lo=0 is used
      phi_lo  <- 0
      Phi_lo  <- 0
    } else {
      tau_lo <- tau[m]
      evalVect <- (tau_lo-x_to_rho)/sig
      phi_lo <- dnorm(evalVect)
      Phi_lo <- pnorm(evalVect)
    }

    if(m == M) {
      tau_hi <- 0 # only used when phi_hi=0 is used
      phi_hi <- 0
      Phi_hi <- 1
    } else {
      tau_hi <- tau[m+1]
      evalVect <- (tau_hi - x_to_rho)/sig
      phi_hi <- dnorm(evalVect)
      Phi_hi <- pnorm(evalVect)
    }

    delta_Phi <- Phi_hi - Phi_lo
    delta_phi <- phi_hi - phi_lo

    # rho
    gradVect[1] <- gradVect[1] + sum(delta_phi/delta_Phi*x_to_rho*log_x/sig)

    # tau_m     [lo]
    if(m > 0) {
      gradVect[m+1] <- gradVect[m+1] + sum(phi_lo/delta_Phi/sig)
    }

    # tau_{m+1} [hi]
    if(m < M) {
      gradVect[m+2] <- gradVect[m+2] - sum(phi_hi/delta_Phi/sig)
    }

    # s
    leadingTerm <- (phi_hi*(tau_hi-x_to_rho)-phi_lo*(tau_lo-x_to_rho))/delta_Phi/sig_sq
    if(hetero) {
      gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm * (1+kap*x))
    } else {
      gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm)
    }

    if(hetero) {
      # kappa
      gradVect[M+3] <- gradVect[M+3] + sum(leadingTerm * s * x)

    }
  }
  return(gradVect)
}

#' @export
fitPowLawOrd <- function(x,v,hetero=F,returnJustParam=T) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  M <- length(unique(v)) - 1
  x_list <- powLawOrdCalc_x_list(x,v)

  rho0   <- 1
  s0     <- mean(x)
  tau0_1 <- s0*qnorm(1/(M+1)) + min(x)
  if(M == 1) {
    tau0 <- tau0_1
  } else {
    tau0_M <- s0*qnorm(M/(M+1)) + max(x)
    dtau0  <- (tau0_M-tau0_1)/(M-1)
    tau0   <- tau0_1 + dtau0*(0:(M-1))
  }

  th_v0 <- c(rho0,tau0,s0)
  if(hetero) {
    kap0   <- 1
    th_v0 <- c(th_v0,kap0)
  }

  optimControl <- list(reltol=1e-12,maxit=100000)
  fit <- optim(th_v0,powLawOrdNegLogLik,gr=powLawOrdGradNegLogLik,method='BFGS',control=optimControl,x_list=x_list,hetero=hetero,hessian=T)

  # By default (returnJustParam=T), return just the optimized parameter vector
  if(returnJustParam) {
    return(fit$par)
  } else {
    return(list(fit=fit,th_v0=th_v0))
  }
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
