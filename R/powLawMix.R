#' @title Power law with mixed ordinal and continuous observations
#'
#' @description \code{powLawMixNegLogLik} calculates the negative log-likelihood. \code{powLawMixGradNegLogLik} calculates the gradient of the negative log-likelihood. \code{extract_th_v} extracts the parameterization for a single ordinal variable from the parameter vector th_y. \code{extract_th_w} extracts the parameterization for a single continuous variable from the parameter vector th_y. \code{simPowLawMix} creats simulated data for a conditionally independent, mixed model
#'
#' @details We assume a mixture of J ordinal variables as described in yada::powLawOrd and K continuous variables as described in yada::powLaw. The calculation of the mean and scale term on the noise, s, is unchanged for the mixed case. The scaling term, kappa, is can be flexibly specified. The ordering of the complete parameter vector th_y is
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
#' kappa     length(unique(modSpec$hetGroups))
#'
#' The extension of the negative log-likelihod calculation from the single variables case ordinal/continuous case to the multi-variable mixed case is straightforward: each variable contributes independently to the sum. Extension of the gradient calculation is similarly straightforward: aside from kappa, which is flexibly specified, the calculation is unchanged. For kappa, a sum across variables is needed.
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
#' @param modSpec Model specification that, among other things, specifies how to unpack th_y
#' @param th_y Parameter vector with ordering th_y = [rho,tau,a,r,b,s,kappa]
#' @param hetero Whether the model is heteroskedastic [Default FALSE]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
powLawMixNegLogLik <- function(th_y,x,Y,modSpec,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- is_hetero(modSpec)

  J <- get_J(modSpec) # number of ordinal variables
  K <- get_K(modSpec) # number of continuous variables

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
    th_v <- extract_th_v(th_y,modSpec,j)
    negLogLik <- negLogLik + powLawOrdNegLogLik(th_v,xj,vj,hetero,transformVar)
  }
 }

  # Add contribution of continuous variables
 if(K > 0) {
  for(k in 1:K) {
    wk <- Y[J+k,]
    indk <- !is.na(x) & !is.na(wk)
    xk <-  x[indk]
    wk <- wk[indk]
    th_w <- extract_th_w(th_y,modSpec,k)
    negLogLik <- negLogLik + powLawNegLogLik(th_w,xk,wk,transformVar)
  }
 }
  return(negLogLik)
}

#' @export
powLawMixGradNegLogLik <- function(th_y,x,Y,modSpec,transformVar=F) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
  # eta_y is the negative log-likelihood
  # For optimization, th_y is the first input
  hetero <- is_hetero(modSpec)
  J <- get_J(modSpec) # number of ordinal variables
  K <- get_K(modSpec) # number of continuous variables

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
    th_v <- extract_th_v(th_y,modSpec,j)
    eta_v <- powLawOrdGradNegLogLik(th_v,xj,vj,hetero,transformVar)
    
    l <- get_var_index('rho',modSpec,j=j)
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[1]
    l <- get_var_index('tau',modSpec,j=j)
    Mj <- modSpec$M[j]
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[2:(1+Mj)]
    l <- get_var_index('s',modSpec,j=j)
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[2+Mj]
    if(hetero) {
      l <- get_var_index('kappa',modSpec,j=j)
      gradNegLogLik[l] <- gradNegLogLik[l] + eta_v[3+Mj]
    }
  }
 
  # Add contribution of continuous variables
  for(k in 1:K) {
    wk <- Y[J+k,]
    indk <- !is.na(x) & !is.na(wk)
    xk <-  x[indk]
    wk <- wk[indk]
    th_w <- extract_th_w(th_y,modSpec,k)
    eta_w <- powLawGradNegLogLik(th_w,xk,wk,transformVar)
    
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kappa]
    l <- get_var_index('a',modSpec,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[1]
    l <- get_var_index('r',modSpec,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[2]
    l <- get_var_index('b',modSpec,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[3]
    l <- get_var_index('s',modSpec,k=k)
    if(l == 1) {
      print(k)
    }
    gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[4]
    if(hetero) {
      l <- get_var_index('kappa',modSpec,k=k)
      gradNegLogLik[l] <- gradNegLogLik[l] + eta_w[5]
    }
  }
  return(gradNegLogLik)
}

#' @export
extract_th_v <- function(th_y,modSpec,j) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,z,kappa]
  # It is assumed that the model is conditionally independent (no z)
  hetero <- is_hetero(modSpec)
  rho <- th_y[get_var_index('rho',modSpec,j=j)]
  tau <- th_y[get_var_index('tau',modSpec,j=j)]
  s   <- th_y[get_var_index('s'  ,modSpec,j=j)]

  if(hetero) {
    kappa <- th_y[get_var_index('kappa',modSpec,j=j)]
  } else {
    kappa <- c()
  }
  return(c(rho,tau,s,kappa))
}

#' @export
extract_th_w <- function(th_y,modSpec,k) {
  # th_y has ordering th_y = [rho,tau,a,r,b,s,z,kappa]
  # It is assumed that the model is conditionally independent (no z)
  hetero <- is_hetero(modSpec)
  a <- th_y[get_var_index('a',modSpec,k=k)]
  r <- th_y[get_var_index('r',modSpec,k=k)]
  b <- th_y[get_var_index('b',modSpec,k=k)]
  s <- th_y[get_var_index('s',modSpec,k=k)]

  if(hetero) {
    kappa <- th_y[get_var_index('kappa',modSpec,k=k)]
  } else {
    kappa <- c()
  }
  return(c(a,r,b,s,kappa))
}

#' @export
simPowLawMix <- function(th_y_list,th_x_list,N,modSpec) {
  if(th_x_list$fitType == 'uniform') {
    x <- runif(N,th_x_list$xmin,th_x_list$xmax)
  } else {
    stop('Only uniform currently supported')
  }

  th_y_vect <- theta_y_list2vect(th_y_list)
  hetero <- is_hetero(modSpec) 
  J <- get_J(modSpec)
  K <- get_K(modSpec)
  Ystar <- matrix(NA,J+K,N)
  for(n in 1:N) {
    for(j in 1:J) {
      th_v <- extract_th_v(th_y_vect,modSpec,j)
      Ystar[j,n] <- rnorm(1,powLawOrd(x[n],th_v),powLawOrdSigma(x[n],th_v,hetero))
    }
    for(k in 1:K) {
      th_w <- extract_th_w(th_y_vect,modSpec,k)
      Ystar[J+k,n] <- rnorm(1,powLaw(x[n],th_w),powLawSigma(x[n],th_w))
    }
  }


  Y <- Ystar
  for(j in 1:J) {
    for(n in 1:N) {
      Y[j,n] <- as.numeric(cut(Ystar[j,n],c(-Inf,th_y_list$tau[[j]],Inf))) - 1 # The ordinal observation
    }
  }
  return(list(x=x,Ystar=Ystar,Y=Y))
}
