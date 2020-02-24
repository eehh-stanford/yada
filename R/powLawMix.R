#' @title Power law with mixed ordinal and continuous observations
#'
#' @description \code{powLawMixNegLogLik} calculates the negative log-likelihood. \code{powLawMixGradNegLogLik} calculates the gradient of the negative log-likelihood. \code{extract_th_v} extracts the parameterization for a single ordinal variable from the parameter vector th_y. \code{extract_th_w} extracts the parameterization for a single continuous variable from the parameter vector th_y.
#'
#' @details We assume a mixture of J ordinal variables as described in yada::powLawOrd and K continuous variables as described in yada::powLaw. The calculation of the mean and scale term on the noise, s, is unchanged for the mixed case. However, the scaling term, kappa, is common across variables. The ordering of the complete parameter vector th_y is
#'
#' \deqn{th_y = [rho,tau,a,r,b,s,kappa]}
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
#' kappa     1
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
#' @param kappa Slope of noise [Optional]
#' @param hp Hyperparameters that, among other things, specify how to unpack th_y
#' @param th_y Parameter vector with ordering th_y = [rho,tau,a,r,b,s,kappa]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLawMixNegLogLik2 <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- is_hetero(hp$paramModel)
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
      kappa <- th_v[3+Mj]
    }
    for(n in 1:N) {
      vnj <- Y[j,n]
      if(!is.na(vnj)) {
        xn   <- x[n]
        gn   <- xn^rho
        sig_n <- s
        if(hetero) {
          sig_n <- sig_n * (1 + kappa*xn)
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
      kappa <- th_w[5]
    }
    for(n in 1:N) {
      wnj <- Y[J+k,n]
      if(!is.na(wnj)) {
        xn   <- x[n]
        hn   <- a*xn^r + b
        sig_n <- s
        if(hetero) {
          sig_n <- sig_n * (1 + kappa*xn)
        }

      negLogLik <- negLogLik - dnorm(wnj,hn,sig_n,log=T)
      }
    }   
  }
  return(negLogLik)
}

#' @export
# hetero is unneeded since hp is given
powLawMixNegLogLik <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- is_hetero(hp$paramModel)
  J <- hp$J # number of ordinal variables
  K <- hp$K # number of continuous variables

  negLogLik <- 0

  # Add contribution of ordinal variables
 if(J > 0) {
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
 }

  # Add contribution of continuous variables
 if(K > 0) {
  for(k in 1:K) {
    wk <- Y[J+k,]
    indk <- !is.na(x) & !is.na(wk)
    xk <-  x[indk]
    wk <- wk[indk]
    th_w <- extract_th_w(th_y,hp,k)
    negLogLik <- negLogLik + powLawNegLogLik(th_w,xk,wk,transformVar)
    
  }
 }
  return(negLogLik)
}

#' @export
# hetero is unneeded since hp is given
powLawMixGradNegLogLik <- function(th_y,x,Y,hp,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- is_hetero(hp$paramModel)
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
      l <- get_var_index('kappa',hp)
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
    eta_w <- powLawGradNegLogLik(th_w,xk,wk,transformVar)
    
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
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
      l <- get_var_index('kappa',hp)
    if(l == 1) {
      print(k)
    }
      gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[5]
    }
  }
  return(gradNegLogLik)
}

#' @export
extract_th_v <- function(th_y,hp,j) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  J <- hp$J
  K <- hp$K
  hetero <- is_hetero(hp$paramModel)
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
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  J <- hp$J
  K <- hp$K
  hetero <- is_hetero(hp$paramModel)
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
simPowLawMix <- function(th_y_list,th_x_list,N,hp) {
  if(th_x_list$fitType == 'uniform') {
    x <- runif(N,th_x_list$xmin,th_x_list$xmax)
  } else {
    stop('Only uniform currently supported')
  }

  th_y_vect <- theta_y_list2vect(th_y_list)
  hetero <- is_hetero(hp$paramModel) 
  J <- hp$J
  K <- hp$K
  Ystar <- matrix(NA,J+K,N)
  for(n in 1:N) {
    for(j in 1:J) {
      th_v <- extract_th_v(th_y_vect,hp,j)
      Ystar[j,n] <- rnorm(1,powLawOrd(x[n],th_v),powLawOrdSigma(x[n],th_v,hetero))
    }
    for(k in 1:K) {
      th_w <- extract_th_w(th_y_vect,hp,k)
      Ystar[J+k,n] <- rnorm(1,powLaw(x[n],th_w),powLawSigma(x[n],th_w))
    }
  }


  Y <- Ystar
  for(j in 1:hp$J) {
    for(n in 1:N) {
      Y[j,n] <- as.numeric(cut(Ystar[j,n],c(-Inf,th_y_list$tau[[j]],Inf))) - 1 # The ordinal observation
    }
  }
  return(list(x=x,Ystar=Ystar,Y=Y))
}
