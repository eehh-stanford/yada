#' @title Transformation functions for theta_y
#'
#' @description \code{powLawOrd}
#'
#' @details Stuff
#' @param th_v Vector of parameters with ordering [rho,tau_1,...,tau_M,s,kap]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
check_model <- function(paramModel) {
  if( !(paramModel %in% known_models()) ) {
    stop(paste('Unrecognized model',paramModel))
  }
}

#' @export
known_models <- function() {
  knownModels <- c('powLawHomo',            # single variable continuous
                   'powLawHetero',          # single variable continuous
                   'powLawOrdHomo',         # single variable ordinal
                   'powLawOrdHetero',       # single variable ordinal
                   'powLawMixUncorrHomo',   # mixed / uncorrelated /   homoskedastic
                   'powLawMixUncorrHetero', # mixed / uncorrelated / heteroskedastic
                   'powLawMixCorrHomo',     # mixed /   correlated /   homoskedastic
                   'powLawMixCorrHetero')   # mixed /   correlated / heteroskedastic
  return(knownModels)
}

#' @export
is_hetero <- function(paramModel) {
  grepl('hetero', tolower(paramModel))
}

#' @export
is_corr <- function(paramModel) {
  # Return the appopriate value explicitly for known model types. Otherwise,
  # through an error.

  if(paramModel == 'powLawHomo') {
    return(F)
  } else if(paramModel == 'powLawHetero') {
    return(F)
  } else if(paramModel == 'powLawOrdHomo') {
    return(F)
  } else if(paramModel == 'powLawOrdHetero') {
    return(F)
  } else if(paramModel == 'powLawMixUncorrHomo') {
    return(F)
  } else if(paramModel == 'powLawMixUncorrHetero') {
    return(F)
  } else if(paramModel == 'powLawMixCorrHomo') {
    return(T)
  } else if(paramModel == 'powLawMixCorrHetero') {
    return(T)
  } else {
    stop(paste('Unrecognized model',paramModel))
  }
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
    indz <- get_var_index('z',hp)
    zbar <- th_y_vect[indz]

    z <- zbar
    N <- (-1 + sqrt(1 + 8*length(z)))/2 # dimension of Sigma
    ind <- (1:N)*(2:(N+1))/2
    z[ind] <- exp(z[ind])
    # build Sigma column by column
    offset <- 0
    U <- matrix(0,N,N)
    for(cc in 1:N) {
      ind_cc <- offset + 1:cc
      U[1:cc,cc] <- z[ind_cc]
      offset <- offset + cc
    }
    th_y_list$Sigma <- t(U) %*% U
  }

  # b is unconstrained. No transformation needed

  # kap should be positive 
  if('kap' %in% names(th_y_list)) {
    th_y_list$kap <- exp(th_y_list$kap)
  }

  return(theta_y_list2vect(th_y_list))
}

#' @export
theta_y_vect2list <- function(th_y_vect,hp) {
  check_model(hp$paramModel)
  hetero <- is_hetero(hp$paramModel)
  corr   <- is_corr  (hp$paramModel)

  if('J' %in% names(hp)) {
    J <- hp$J
  } else {
    J <- 0
  }

   if('K' %in% names(hp)) {
    K <- hp$K
  } else {
    K <- 0
  }
 
  th_y_list <- list(paramModel=hp$paramModel)
  if(!corr) {
    s   <- rep(NA,J+K) # add s later to keep ordering of list variables
  }
  # For correlated models, Sigma is extracted below

  if(J > 0) {
    th_y_list$rho <- rep(NA,J)
    th_y_list$tau <- list()
    for(j in 1:J) {
      th_y_list$rho[ j ] <- th_y_vect[get_var_index('rho',hp,j=j)]
      th_y_list$tau[[j]] <- th_y_vect[get_var_index('tau',hp,j=j)]
      if(!corr) {
        s[ j ] <- th_y_vect[get_var_index('s',  hp,j=j)]
      }
    }
  }

  if(K > 0) {
    th_y_list$a <- rep(NA,K)
    th_y_list$r <- rep(NA,K)
    th_y_list$b <- rep(NA,K)
    for(k in 1:K) {
      th_y_list$a[k]   <- th_y_vect[get_var_index('a',hp,k=k)]
      th_y_list$r[k]   <- th_y_vect[get_var_index('r',hp,k=k)]
      th_y_list$b[k]   <- th_y_vect[get_var_index('b',hp,k=k)]
      if(!corr) {
        s[J+k] <- th_y_vect[get_var_index('s',hp,k=k)]
      }
    }
  }

  if(!corr) {
    th_y_list$s <- s
  } else {
    z <- th_y_vect[get_var_index('z',hp)]
    U <- matrix(0,nrow=J+K,ncol=J+K)
    U[upper.tri(U,diag=T)] <- z
    th_y_list$Sigma <- t(U) %*% U
  }

  if(hetero) {
    th_y_list$kap <- th_y_vect[get_var_index('kappa',hp)]
  }

  return(th_y_list)
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
    # Either s or Sigma is a field
    if('s' %in% names(th_y_list)) {
      s <- th_y_list$s
    } else {
      s <- sqrt(diag(th_y_list$Sigma))
    }
    th_y_vect <- c(th_y_vect,s)

  } else {
    # Cholesky decomposition (upper triangular)
    U <- chol(th_y_list$Sigma)

    # Get upper diagonal elements of U
    # Unwraps by column: c(reduced_col1,reduced_col2,reduced_col3,...)
    th_y_vect <- c(th_y_vect,U[upper.tri(U,diag=T)])
  }

  if(hetero) {
    th_y_vect <- c(th_y_vect,th_y_list$kap)
  }

  return(th_y_vect)
}

#' @export
get_var_index <- function(varName,hp,j=NA,k=NA) {
  # This is a helper function to return a variable's index given the variable
  # index j and variable name. This task is centralized here in large part to
  # improve code readability. For tau, a vector is returned
  #
  # th_y has ordering th_y = [rho,tau,a,r,b,s,kap]
  check_model(hp$paramModel)
  hetero <- is_hetero(hp$paramModel)
  corr   <- is_corr  (hp$paramModel)

  if('J' %in% names(hp)) {
    J <- hp$J
  } else {
    J <- 0
  }

   if('K' %in% names(hp)) {
    K <- hp$K
  } else {
    K <- 0
  }
 
  if(varName == 'kappa') {
    if(!hetero) {
      stop('kappa requested but model is not heteroskedastic')
    }

    if(!corr) {
      offset <- 2*J + sum(hp$M) + 4*K
    } else {
      offset <- 2*J + sum(hp$M) + 4*K + choose(J+K,2)
    }
    return(offset+1)
  }

  if(varName == 'z') {
    if(!corr) {
      stop('z requested but model is not correlated')
    }

    offset <- 1*J + sum(hp$M) + 3*K
    return((offset+1):(offset+J+K+choose(J+K,2)))
  }

  if(!is.na(j) && !(is.na(k))) {
    stop('Both j and k are specified')
  }

  if(is.na(j) && (is.na(k))) {
    stop('Neither j nor k is specified')
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


