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
      s   <- exp(th_v[numParam-1])
      kap <- exp(th_v[numParam])
    }
    sig <- s*(1 + kap*x)
  } else {
    if(transformVar) { 
      s   <- exp(th_v[numParam-1])
    } else {
      s   <- th_v[numParam-1]
    }
    sig <- s
  }

  return(sig)
}

#' @export
powLawOrdNegLogLik <- function(th_v,x_list,hetero=F,transformVar=F,hp=NA) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # eta_v is the negative log-likelihood
  # For optimization, make th_v the first input

  if(transformVar) {
    # hp is only needed if transformVar is T, and this dependency could be avoided
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
fitPowLawOrd <- function(x,v,hetero=F,returnJustParam=T,transformVar=F) {
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

  if(transformVar) {
    rho0 <- log(rho0)
    s0   <- log(s0)
  }
  th_v0 <- c(rho0,tau0,s0)
  if(hetero) {
    kap0   <- 1
    if(transformVar) {
      kap0 <- log(kap0)
    }
    th_v0 <- c(th_v0,kap0)
  }

  optimControl <- list(reltol=1e-12,maxit=100000)
  #optimControl <- list(maxit=100000,trace=6)
  #fit <- optim(th_v0,powLawOrdNegLogLik,control=optimControl,x_list=x_list,hetero=hetero,hessian=T,transformVar=transformVar)
  fit <- optim(th_v0,powLawOrdNegLogLik,gr=powLawOrdGradNegLogLik,method='BFGS',control=optimControl,x_list=x_list,hetero=hetero,hessian=T,transformVar=transformVar)
  #lowerBounds <- c(-Inf,rep(-Inf,M),0)
  #upperBounds <- c( Inf,rep( Inf,M),Inf)
  #if(hetero) {
  #  lowerBounds <- c(lowerBounds,0)
  #  upperBounds <- c(lowerBounds,Inf)
  #}
  #fit <- optim(th_v0,powLawOrdNegLogLik,gr=powLawOrdGradNegLogLik,method='L-BFGS-B',control=optimControl,x_list=x_list,hetero=hetero,hessian=T,lower=lowerBounds,upper=upperBounds)

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
theta_y_constr2unconstr <- function(th_y_vect,hp) {
  check_model(hp$paramModel)

  # For code clarity, convert to a list representation
  th_y_list <- theta_y_vect2list(th_y_vect,hp)

  # rho should be positive 
  if('rho' %in% names(th_y_list)) {
    th_y_list$rho <- log(th_y_list$rho)
  }

  # tau_1 is unconstrained. Successive differences should be positive
  if('tau' %in% names(th_y_list)) {
    for(j in 1:length(th_y_list$tau)) {
      tau_j <- th_y_list$tau[[j]]
      M_j   <- length(tau_j)
      th_y_list$tau[[j]] <- tau_j[1]
      if(M_j > 1) {
        th_y_list$tau[[j]] <- c(th_y_list$tau[[j]],log(tau_j[2:M_j] - tau_j[1:(M_j-1)]))
      }
    }
  }

  # a should be positive 
  if('a' %in% names(th_y_list)) {
    th_y_list$a <- log(th_y_list$a)
  }

  # r should be positive 
  if('r' %in% names(th_y_list)) {
    th_y_list$r <- log(th_y_list$r)
  }

  # s should be positive 
  if('s' %in% names(th_y_list)) {
    th_y_list$s <- log(th_y_list$s)
  }

  # handle Sigma a little differently
  # Only the diagonal elements of U need to be positive
  if('Sigma' %in% names(th_y_list)) {
    #th_y_list$s <- log(th_y_list$s)
    indz <- get_var_index('z',hp)
    z    <- th_y_vect[indz]
    zbar <- z
    U <- matrix(0,nrow=hp$J+hp$K,ncol=hp$J+hp$K)
    U[upper.tri(U,diag=T)] <- z
    #th_y_list$Sigma <- t(U) %*% U # This is effectively ignored
    N <- hp$J+hp$K
    #ind <- (1:N-1)*N - (1:N-1)*(1:N-2)/2 + 1
    ind <- (1:N)*(2:(N+1))/2
    zbar[ind] <- log(zbar[ind])
    th_y_vect[indz] <- zbar
  }

  # b is unconstrained. No transformation needed

  # kap should be positive 
  if('kap' %in% names(th_y_list)) {
    th_y_list$kap <- log(th_y_list$kap)
  }

  th_y_vect <- theta_y_list2vect(th_y_list)

  if('Sigma' %in% names(th_y_list)) {
    th_y_vect[indz] <- zbar
  }
  
  return(th_y_vect)
}

#' @export
theta_y_unconstr2constr <- function(th_y_vect,hp) {
  check_model(hp$paramModel)

  # For code clarity, convert to a list representation
  th_y_list <- theta_y_vect2list(th_y_vect,hp)

  # rho should be positive 
  if('rho' %in% names(th_y_list)) {
    th_y_list$rho <- exp(th_y_list$rho)
  }

  # tau_1 is unconstrained. Successive differences should be positive
  if('tau' %in% names(th_y_list)) {
    for(j in 1:length(th_y_list$tau)) {
      tau_j <- th_y_list$tau[[j]]
      M_j   <- length(tau_j)
      if(M_j == 1) {
        th_y_list$tau[[j]] <- tau_j
      } else {
        th_y_list$tau[[j]] <- tau_j[1] + c(0,cumsum(exp(tau_j[2:M_j])))
      }
    }
  }

  # a should be positive 
  if('a' %in% names(th_y_list)) {
    th_y_list$a <- exp(th_y_list$a)
  }

  # r should be positive 
  if('r' %in% names(th_y_list)) {
    th_y_list$r <- exp(th_y_list$r)
  }

  # s should be positive 
  if('s' %in% names(th_y_list)) {
    th_y_list$s <- exp(th_y_list$s)
  }

  if('Sigma' %in% names(th_y_list)) {
    stop('Implement for corr')
    #Sigma <- th_y_list$Sigma
    #Sigma[1 + (nrow(Sigma)+1)*(0:(nrow(Sigma)-1))] <- exp(diag(Sigma))
    #th_y_list$Sigma <- Sigma
  }

  # b is unconstrained. No transformation needed

  # kap should be positive 
  if('kap' %in% names(th_y_list)) {
    th_y_list$kap <- exp(th_y_list$kap)
  }

  return(theta_y_list2vect(th_y_list))
}

#' @export
theta_y_list2vect <- function(th_y_list) {
  check_model(th_y_list$paramModel)
  hetero <- is_hetero(th_y_list$paramModel)
  corr   <- is_corr  (th_y_list$paramModel)

  if('rho' %in% names(th_y_list)) {
    J <- length(th_y_list$rho)
  } else {
    J <- 0
  }

  if('r' %in% names(th_y_list)) {
    K <- length(th_y_list$r)
  } else {
    K <- 0
  }

  th_y_vect <- c()
  if(J > 0) {
    th_y_vect <- c(th_y_vect,th_y_list$rho)
    th_y_vect <- c(th_y_vect,unlist(th_y_list$tau))
  }

  if(K > 0) {
    th_y_vect <- c(th_y_vect,th_y_list$a)
    th_y_vect <- c(th_y_vect,th_y_list$r)
    th_y_vect <- c(th_y_vect,th_y_list$b)
  }

  if(!corr) {
    th_y_vect <- c(th_y_vect,th_y_list$s)
  } else {
    # Cholesky decomposition (upper triangular)
    U <- chol(th_y_list$Sigma)

    # Get upper diagonal elements of U
    # Unwraps by column: c(col1,reduced_col2,reduced_col3,...)
    th_y_vect <- c(th_y_vect,th_y_list$s,U[upper.tri(U,diag=T)])
  }

  if(hetero) {
    th_y_vect <- c(th_y_vect,th_y_list$kap)
  }

  return(th_y_vect)
}

#' @export
check_model <- function(paramModel) {
  if( !(paramModel %in% known_models()) ) {
    stop(paste('Unrecognized model',paramModel))
  }
}
