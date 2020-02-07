#' @title Power law with ordinal observations
#'
#' @description \code{powLawMix} calculates the mean (h). \code{powLawMixSigma} calculates the noise (sigma, or sig for short). \code{powLawMixNegLogLik} calculates the negative log-likelihood. \code{powLawMixGradNegLogLik} calculates the gradient of the negative log-likelihood. \code{fitPowLawMix} returns the maximum likelihood fit.
#'
#' @details We assume a mixture of J ordinal variables as described in yada::powLawOrd and K continuous variables as described in yada::powLaw. The calculation of the mean and scale term on the noise, s, is unchanged for the mixed case. However, the scaling term, kappa (kap for short, is common across variables. The ordering of the complete parameter vector th_y is
#'
#' \deqn{th_y = [rho,tau,a,r,b,s,kap]}
#'
#' The length of each vector in th_y is
#'
#' Variable  Length
#' rho       J
#' tau       M1 + M2 + ... + Mj + ... + MJ, where Mj+1 is the number of ordinal
#'           categories for the j-th ordinal variable
#' a         K
#' r         K
#' b         K
#' s         J + K
#' kap       1
#'
#' The extension of the negative log-likelihod calculation from the single variables case ordinal/continuous case to the multi-variable mixed case is straightforward: each variable contributes independently to the sum. Extension of the gradient calculation is similarly straightforward: aside from kappa, which is shared across variables, the calculation is unchanged. For kappa, a sum across variables is needed.
#' 
#' @param x Vector of independent variable observations
#' @param vstar Vector of latent dependent variable observations (for v^*)
#' @param v Vector of dependent variable observations (ordinal)
#' @param w Vector of dependent variable observations (continuous)
#' @param rho Scaling exponent (ordinal)
#' @param a Multiplicative coefficient (continuous)
#' @param r Scaling exponent (continuous)
#' @param b Offset (continuous)
#' @param s Baseline noise
#' @param kap Slope of noise (short for kappa) [Optional]
#' @param hp Hyperparameters that, among other things, specify how to unpack th_y
#' @param th_y Parameter vector with ordering th_y = [rho,tau,a,r,b,s,kap]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
# hetero is unneeded since hp is given
powLawMixNegLogLik <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- !grepl('homosk',hp$paramModel)
  J <- hp$J # number of ordinal variables
  K <- hp$K # number of continuous variables

  negLogLik <- 0

  # Add contribution of ordinal variables
  for(j in 1:J) {
    # Extracting xList prior and making it an input to powLawMixNegLogLik
    # would speed up computation. However, extracting it here likely makes the
    # code easier to understand.
    vj <- Y[j,]
    indj <- !is.na(x) & !is.na(vj)
    xj <-  x[indj]
    vj <- vj[indj]
    x_list <- powLawOrdCalc_x_list(xj,vj)
    th_v <- extract_th_v(th_y,hp,j)
    negLogLik <- negLogLik + powLawOrdNegLogLik(th_v,x_list,hetero,transformVar)
  }

  # Add contribution of continuous variables
  for(k in 1:K) {
    wk <- Y[J+k,]
    indk <- !is.na(x) & !is.na(wk)
    xk <-  x[indk]
    wk <- wk[indk]
    th_w <- extract_th_w(th_y,hp,k)
    negLogLik <- negLogLik + powLawNegLogLik(th_w,xk,wk,hetero)
    
  }

  return(negLogLik)
}

#' @export
powLawMixNegLogLik2 <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- !grepl('homosk',hp$paramModel)
  J <- hp$J # number of ordinal variables
  K <- hp$K # number of continuous variables

  negLogLik <- 0

  N <- length(x)
  for(j in 1:J) {
    th_v <- extract_th_v(th_y,hp,j)
    Mj <- hp$M[j]
    rho <- th_v[1]
    tau <- th_v[2:(1+Mj)]
    s   <- th_v[2+Mj]
    if(hetero) {
      kap <- th_v[3+Mj]
    }
    for(n in 1:N) {
      vnj <- Y[j,n]
      if(!is.na(vnj)) {
        xn   <- x[n]
        gn   <- xn^rho
        sig_n <- s
        if(hetero) {
          sig_n <- sig_n * (1 + kap*xn)
        }

        if(vnj == 0) {
          Phi_lo <- 0
        } else {
          tau_lo <- tau[vnj]
          Phi_lo <- pnorm( (tau_lo - gn)/sig_n )
        }

        if(vnj == Mj) {
          Phi_hi <- 1
        } else {
          tau_hi <- tau[vnj+1]
          Phi_hi <- pnorm( (tau_hi - gn)/sig_n )
        }
      negLogLik <- negLogLik - log(Phi_hi - Phi_lo)
      }
    }
  }

  for(k in 1:K) {
    th_w <- extract_th_w(th_y,hp,k)
    a <- th_w[1]
    r <- th_w[2]
    b <- th_w[3]
    s <- th_w[4]

    if(hetero) {
      kap <- th_w[5]
    }
    for(n in 1:N) {
      wnj <- Y[J+k,n]
      if(!is.na(wnj)) {
        xn   <- x[n]
        hn   <- a*xn^r + b
        sig_n <- s
        if(hetero) {
          sig_n <- sig_n * (1 + kap*xn)
        }

      negLogLik <- negLogLik - dnorm(wnj,hn,sig_n,log=T)
      }
    }   
  }
  return(negLogLik)
}

#' @export
# hetero is unneeded since hp is given
powLawMixGradNegLogLik <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- !grepl('homosk',hp$paramModel)
  J <- hp$J # number of ordinal variables
  K <- hp$K # number of continuous variables

  gradNegLogLik <- rep(0,length(th_y))

  # Add contribution of ordinal variables
  for(j in 1:J) {
    # Extracting xList prior and making it an input to powLawMixNegLogLik
    # would speed up computation. However, extracting it here likely makes the
    # code easier to understand.
    vj <- Y[j,]
    indj <- !is.na(x) & !is.na(vj)
    xj <-  x[indj]
    vj <- vj[indj]
    x_list <- powLawOrdCalc_x_list(xj,vj)
    th_v <- extract_th_v(th_y,hp,j)
    eta_v <- powLawOrdGradNegLogLik(th_v,x_list,hetero,transformVar)
    
    l <- get_var_index('rho',hp,j=j)
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[1]
    l <- get_var_index('tau',hp,j=j)
    Mj <- hp$M[j]
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[2:(1+Mj)]
    l <- get_var_index('s',hp,j=j)
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[2+Mj]
    if(hetero) {
      l <- get_var_index('kap',hp)
      gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[3+Mj]
    }
  }
 
  # Add contribution of continuous variables
  for(k in 1:K) {
    wk <- Y[J+k,]
    indk <- !is.na(x) & !is.na(wk)
    xk <-  x[indk]
    wk <- wk[indk]
    th_w <- extract_th_w(th_y,hp,k)
    eta_w <- powLawGradNegLogLik(th_w,xk,wk,hetero)
    
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
    l <- get_var_index('a',hp,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[1]
    l <- get_var_index('r',hp,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[2]
    l <- get_var_index('b',hp,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[3]
    l <- get_var_index('s',hp,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[4]
    if(hetero) {
      l <- get_var_index('kap',hp)
    if(l == 1) {
      print(k)
    }
      gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[5]
    }
  }
  return(gradNegLogLik)
}

#' @export
get_var_index <- function(varName,hp,j=NA,k=NA) {
  # This is a helper function to return a variable's index given the variable
  # index j and variable name. This task is centralized here in large part to
  # improve code readability. For tau, a vector is returned
  #
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]

  J <- hp$J
  K <- hp$K
  if(varName == 'kap') {
    # No error is thrown if the model is homoskedastic
    offset <- 2*J + sum(hp$M) + 4*K
    return(offset+1)
  }

  if(!is.na(j) && !(is.na(k))) {
    stop('Both j and k are specified')
  }

  if(is.na(j) && (is.na(k))) {
    stop('Neither j nor is specified')
  }


  if(!is.na(j)) {
    if(!(varName %in% c('rho','tau','s')) ) {
      stop('Unsupported variable for j being specified')
    }
    if(j < 1 || J < j) {
      stop('j is not between 1 and J')
    }
  }

  if(!is.na(k)) {
    if(!(varName %in% c('a','r','b','s')) ) {
      stop('Unsupported variable for k being specified')
    }
    if(k < 1 || K < k) {
      stop('k is not between 1 and K')
    }
  }

  if(varName == 'rho') {
    return(j)
  } else if(varName == 'tau') {
    if(j == 1) {
      offset <- J
    } else {
      offset <- J + sum(hp$M[1:(j-1)])
    }
    return((offset + 1):(offset+hp$M[j]))
  } else if(varName == 'a') {
    offset <- J + sum(hp$M)
    return(offset + k)
  } else if(varName == 'r') {
    offset <- J + sum(hp$M) + K
    return(offset + k)
  } else if(varName == 'b') {
    offset <- J + sum(hp$M) + 2*K
    return(offset + k)
  } else if(varName == 's') {
    if(!is.na(j)) {
      offset <- J + sum(hp$M) + 3*K
      return(offset + j)
    } else {
      offset <- 2*J + sum(hp$M) + 3*K
      return(offset + k)
    }
  } else {
    stop(paste('Unrecognized variable',varName))
  }
}

#' @export
extract_th_v <- function(th_y,hp,j) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  J <- hp$J
  K <- hp$K
  hetero <- !grepl('homosk',hp$paramModel)
  th_v <- th_y[j] # add rho

  if(j == 1) {
    lastIndex <- J
  } else {
    lastIndex <- J + sum(hp$M[1:(j-1)])
  }
  th_v <- c(th_v,th_y[(lastIndex+1):(lastIndex+hp$M[j])]) # add tau

  lastIndex <- J + sum(hp$M) + 3*K

  th_v <- c(th_v,th_y[lastIndex+j]) # add s

  if(hetero) {
    th_v <- c(th_v,th_y[length(th_y)]) # add kappa
  }
  return(th_v)
}

#' @export
extract_th_w <- function(th_y,hp,k) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  J <- hp$J
  K <- hp$K
  hetero <- !grepl('homosk',hp$paramModel)
  lastIndex <- J + sum(hp$M)
  th_w <- th_y[lastIndex+k] # add a
  lastIndex <- J + sum(hp$M) + K
  th_w <- c(th_w,th_y[lastIndex+k]) # add r
  lastIndex <- J + sum(hp$M) + 2*K
  th_w <- c(th_w,th_y[lastIndex+k]) # add b
  lastIndex <- 2*J + sum(hp$M) + 3*K
  th_w <- c(th_w,th_y[lastIndex+k]) # add s

  if(hetero) {
    th_w <- c(th_w,th_y[length(th_y)]) # add kappa
  }
  return(th_w)
}

#' @export
powLawOrdGradNegLogLik <- function(th_v,x_list,hetero=F,transformVar=T) {
  # th_v has ordering [rho,tau_1,...tau_2,s,kap]
  # eta_v is the negative log-likelihood
  # For optimization, th_v is the first input

  M <- length(x_list) - 1 # number of ordinal categories is M+1
  rho <- th_v[1]         # rho
  tau <- th_v[2:(M+1)]   # tau_1 ... tau_M
  s   <- th_v[M+2]       # s
  if(transformVar) {
    rho <- exp(rho)
    s   <- exp(s)
  }

  if(hetero) {
    kap <- th_v[M+3]     # kappa
    if(transformVar) {
      kap <- exp(kap)
    }
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
    if(!transformVar) {
      gradVect[1] <- gradVect[1] + sum(delta_phi/delta_Phi*x_to_rho*log_x/sig)
    } else {
      gradVect[1] <- gradVect[1] + sum(delta_phi/delta_Phi*x_to_rho*log_x/sig)*rho
    }

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
      if(!transformVar) {
        gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm * (1+kap*x))
      } else {
        gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm * (1+kap*x))*s
      }
    } else {
      if(!transformVar) {
        gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm)
      } else {
        gradVect[M+2] <- gradVect[M+2] + sum(leadingTerm)*s
      }
    }

    if(hetero) {
      # kappa
      if(!transformVar) {
        gradVect[M+3] <- gradVect[M+3] + sum(leadingTerm * s * x)
      } else {
        gradVect[M+3] <- gradVect[M+3] + sum(leadingTerm * s * x)*kap
      }

    }
  }
  return(gradVect)
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
