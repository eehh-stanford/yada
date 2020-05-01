#' @title Transformation functions for theta_y
#'
#' @description \code{powLawOrd}
#'
#' @details Stuff
#' @param th_v Vector of parameters with ordering [rho,tau_1,...,tau_M,s,kappa]
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
check_model <- function(modSpec) {
  if(tolower(modSpec$meanSpec) != 'powlaw') {
    stop(paste('Unrecognized specification for the mean,',modSpec$meanSpec))
  }

  # Reject single variable models specified as conditionally independent
  J <- get_J(modSpec)
  K <- get_K(modSpec)
  if(J + K == 1) {
    if('cdepSpec' %in% names(modSpec)) {
      if(tolower(modSpec$cdepSpec) == 'dep') {
        stop('Model is single-variable, but conditional independence is specified')
      }
    }
  }

  # Reject models specified as conditionally dependent for which cdepGroups is
  # mis-specified.
  if('cdepSpec' %in% names(modSpec)) {
    if(tolower(modSpec$cdepSpec) == 'dep') {
      if( !('cdepGroups' %in% names(modSpec)) ) {
        stop('Model is conditionally dependent, but cdepGroups not given')
      }
      N <- length(modSpec$cdepGroups)
      if(N != J + K) {
        stop(paste('Length of cdepGroups =',N,'but should be J+K=',J+K))
      }
      if(all(is.na(modSpec$cdepGroups))) {
        stop(paste('cdepGroups is all NA'))
      }
   
      Gz <- max(modSpec$cdepGroups,na.rm=T)
      uniqueVal <- sort(unique(modSpec$cdepGroups[!is.na(modSpec$cdepGroups)]))
      if(!all(uniqueVal == 1:Gz)) {
        stop('Unique values of cdepGroups is not 1:Gz')
      }
    }
  }

  # Reject models specified as heteroskedatic for which hetGroups is
  # mis-specified.
  if('hetSpec' %in% names(modSpec)) {
    if(tolower(modSpec$hetSpec) == 'linearsd') { # At least for now, this is the only heteroskedastic parameterization
      if( !('hetGroups' %in% names(modSpec)) ) {
        stop('Model is heteroskedastic, but hetGroups not given')
      }
      N <- length(modSpec$hetGroups)
      if(N != J + K) {
        stop(paste('Length of hetGroups =',N,'but should be J+K=',J+K))
      }
      if(all(is.na(modSpec$hetGroups))) {
        stop(paste('hetGroups is all NA'))
      }
   
      Gkappa <- max(modSpec$hetGroups,na.rm=T)
      uniqueVal <- sort(unique(modSpec$hetGroups[!is.na(modSpec$hetGroups)]))
      if(!all(uniqueVal == 1:Gkappa)) {
        stop('Unique values of hetGroups is not 1:Gkappa')
      }
    }
  }
}

#' @export
get_J <- function(modSpec) {
  # A helper function to get J, which is assumed 0 if J is not a field in
  # modSpec.
  if('J' %in% names(modSpec)) {
    return(modSpec$J)
  } else {
    return(0)
  }
}

#' @export
get_K <- function(modSpec) {
  # A helper function to get K, which is assumed 0 if K is not a field in
  # modSpec.
  if('K' %in% names(modSpec)) {
    return(modSpec$K)
  } else {
    return(0)
  }
}

#' @export
get_Gkappa <- function(modSpec) {
  # A helper function to get Gkappa, which is assumed 0 if hetGroups is not a
  # field in modspec
  if('hetGroups' %in% names(modSpec)) {
    return(max(modSpec$hetGroups,na.rm=T))
  } else {
    return(0)
  }
}

#' @export
get_Gz <- function(modSpec) {
  # A helper function to get Gz, which is assumed 0 if cdepGroups is not a
  # field in modspec
  if('cdepGroups' %in% names(modSpec)) {
    return(max(modSpec$cdepGroups,na.rm=T))
  } else {
    return(0)
  }
}

#' @export
is_hetero <- function(modSpec) {
  # The model is homoskedastic if (1) modSpec$hetSpec is 'none' or if
  # (2) modSpec$hetSpec is 'linearSd', but hetGroups is all NA.
  if(tolower(modSpec$hetSpec) == 'none') {
    return(F)
  } else if(tolower(modSpec$hetSpec) == 'linearsd') {
    return(!all(is.na(modSpec$hetGroups)))
  } else {
    stop(paste('Unsupported hetSpec,',modSpec$hetSpec))
  }
}


#' @export
is_cdep <- function(modSpec) {
  # For the special case of a one-variable model, False is returned. Aside from
  # this, the model is conditionally depenendent if (a) modSpec$cdepSpec is
  # 'dep' or (b) modSpec$cdepSpec is 'indep', but cdepGroups is all NA.

  J <- get_J(modSpec)
  K <- get_K(modSpec)
  if(J + K == 1) {
    return(F)
  }

  if(tolower(modSpec$cdepSpec) == 'indep') {
    return(F)
  } else if(tolower(modSpec$cdepSpec) == 'dep') {
    return(!all(is.na(modSpec$cdepGroups)))
  } else {
    stop(paste('Unsupported cdepSpec,',modSpec$cdepSpec))
  }
}

#' @export
get_z_length <- function(modSpec) {
  if( !('cdepSpec' %in% names(modSpec)) ) {
    return(0)
  }
  if(tolower(modSpec$cdepSpec == 'indep')) {
    return(0)
  }
  groupSizes <- as.vector(table(modSpec$cdepGroups))
  numGroups  <- length(as.vector(table(modSpec$cdepGroups)))
  return(sum(groupSizes > 1) + choose(numGroups,2))
}

#' @export
get_non_singleton_groups <- function(groups) {
  return(sort(as.numeric(names(counts))[counts > 1]))
}

#get_num_nonsing_groups <- function(groups) {
#  return(sum(as.vector(table(modSpec$cdepGroups)) > 1))
#}

#get_num_groups <- function(groups) {
#  return(length(unique(groups)))
#}

#' @export
get_var_index <- function(varName,modSpec,j=NA,k=NA,i1=NA,i2=NA) {
  # This is a helper function to return a variable's index given the variable
  # index j or k and variable name. This task is centralized here in large part
  # to improve code readability and make test easier. For tau, a vector is
  # returned.
  #
  # th_y has ordering th_y = [rho,tau,a,r,b,s,z,kappa]


  check_model(modSpec)
  hetero <- is_hetero(modSpec)
  cdep   <- is_cdep  (modSpec)

  J <- get_J(modSpec)
  K <- get_K(modSpec)

  if(varName == 'z') {
    if(!cdep) {
      stop('z requested but model is not conditionally dependent')
    }

    #         rho   M                a/r/b  s
    offset <- 1*J + sum(modSpec$M) + 3*K  + J+K

    if(is.na(i1) && is.na(i2)) {
      # Then return all the indices
      return(offset + (1:get_z_length(modSpec)))
    }

    if(i1 == i2) {
      stop('i1 should not equal i2')
    }

    # If i1 and i2 are members of the same group, this is an intragroup
    # correlation that is stored in the beginning of z.
    g1 <- modSpec$cdepGroups[i1]
    g2 <- modSpec$cdepGroups[i2]
    if(g1 == g2) {
      return(offset + g1)
    }

    # If i1 and i2 are members of different groups, this is an intergroup
    # correlation that is stored after the intragroup correlations.

    offset <- offset + sum(as.vector(table(modSpec$cdepGroups)) > 1)
    numGroups <- length(unique(modSpec$cdepGroups))
    if(g1 < g2) {
      index <- elemToIndex(c(g1-1,g2-1),numGroups) + 1
    } else {
      index <- elemToIndex(c(g2-1,g1-1),numGroups) + 1
    }
    return(offset+index)
  }

  if(varName == 'kappa') {
    if(!hetero) {
      stop('kappa requested but model is not heteroskedastic')
    }

    # If modSpec$M is not in M, sum(modSpec$M) evaluates to 0
    #         rho   M                a/r/b  s     z
    offset <- 1*J + sum(modSpec$M) + 3*K  + J+K + get_z_length(modSpec)

    if(is.na(j) && is.na(k)) {
      # Then return all the indices
      return(offset + (1:length(unique(modSpec$hetGroups))))
    }

    if(!is.na(j)) {
      return(offset+modSpec$hetGroups[j])
    } else {
      return(offset+modSpec$hetGroups[J+k])
    }
  }

  if(!is.na(j) && !(is.na(k))) {
    stop('Both j and k are specified')
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

  # th_y has ordering th_y = [rho,tau,a,r,b,kappa,s,z]
  if(varName == 'rho') {
    if(is.na(j)) {
      return(1:J)
    } else {
      return(j)
    }
  } else if(varName == 'tau') {
    if(is.na(j)) {
      stop('j must be specified for variable tau')
    }
    if(j == 1) {
      offset <- J
    } else {
      offset <- J + sum(modSpec$M[1:(j-1)])
    }
    return((offset + 1):(offset+modSpec$M[j]))
  } else if(varName == 'a') {
    #         rho   M
    offset <- 1*J + sum(modSpec$M)
    if(is.na(k)) {
      return(offset + 1:K)
    } else {
      return(offset + k)
    }
  } else if(varName == 'r') {
    #         rho   M                a
    offset <- 1*J + sum(modSpec$M) + 1*K
    if(is.na(k)) {
      return(offset + 1:K)
    } else {
      return(offset + k)
    }
  } else if(varName == 'b') {
    #         rho   M                a/r
    offset <- 1*J + sum(modSpec$M) + 2*K
    if(is.na(k)) {
      return(offset + 1:K)
    } else {
      return(offset + k)
    }
  } else if(varName == 's') {
    #         rho   M                a/r/b
    offset <- 1*J + sum(modSpec$M) + 3*K
    if(!is.na(j)) {
      return(offset + j)
    } else if(!is.na(k)) {
      return(offset + J + k)
    } else {
      return(offset + 1:(J+K))
    }
  } else {
    stop(paste('Unrecognized variable',varName))
  }
}

#' @export
theta_y_vect2list <- function(th_y_vect,modSpec) {
  check_model(modSpec)
  hetero <- is_hetero(modSpec)
  cdep   <- is_cdep  (modSpec)

  J <- get_J(modSpec)
  K <- get_K(modSpec)
 
  
  # Create the output list, adding the model specification to it
  th_y_list <- list(modSpec=modSpec)

  # Initialize the vector of baseline standard deviations
  s   <- rep(NA,J+K) # add s later to keep ordering of list variables

  if(J > 0) {
    th_y_list$rho <- rep(NA,J)
    th_y_list$tau <- list()
    for(j in 1:J) {
      th_y_list$rho[ j ] <- th_y_vect[get_var_index('rho',modSpec,j=j)]
      th_y_list$tau[[j]] <- th_y_vect[get_var_index('tau',modSpec,j=j)]
      s[j] <- th_y_vect[get_var_index('s',modSpec,j=j)]
    }
  }

  if(K > 0) {
    th_y_list$a <- rep(NA,K)
    th_y_list$r <- rep(NA,K)
    th_y_list$b <- rep(NA,K)
    for(k in 1:K) {
      th_y_list$a[k]   <- th_y_vect[get_var_index('a',modSpec,k=k)]
      th_y_list$r[k]   <- th_y_vect[get_var_index('r',modSpec,k=k)]
      th_y_list$b[k]   <- th_y_vect[get_var_index('b',modSpec,k=k)]
      s[J+k] <- th_y_vect[get_var_index('s',modSpec,k=k)]
    }
  }

  Sigma <- diag(s^2)

  # For a conditionally dependent model, add the off-diagonal terms to Sigma
  if(cdep) {
#    z <- th_y_vect[get_var_index('z',modSpec)]
    for(i1 in 1:(J+K-1)) {
      g1 <- modSpec$cdepGroups[i1]
      for(i2 in (i1+1):(J+K)) {
        g2 <- modSpec$cdepGroups[i2]
        z <- th_y_vect[get_var_index('z',modSpec,i1=i1,i2=i2)]
        Sigma[i1,i2] <- s[i1]*s[i2]*z
        Sigma[i2,i1] <- Sigma[i1,i2]
      }
    }
  }
  th_y_list$Sigma <- Sigma

  # For a heteroskedastic model, build the full-length kappa (length J+K)
  if(hetero) {
    th_y_list$kappa <- rep(NA,J+K)
    for(j in 1:J) {
      th_y_list$kappa[j] <- th_y_vect[get_var_index('kappa',modSpec,j=j)]
    }
    for(k in 1:K) {
      th_y_list$kappa[J+k] <- th_y_vect[get_var_index('kappa',modSpec,k=k)]
    }
  }

  return(th_y_list)
}

#' @export
theta_y_list2vect <- function(th_y_list) {
  modSpec <- th_y_list$modSpec # Extract to improve code readability
  check_model(modSpec)
  # TODO: add a check here of the parameterization consistency
  hetero <- is_hetero(modSpec)
  cdep   <- is_cdep  (modSpec)

  J <- get_J(modSpec)
  K <- get_K(modSpec)
 
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

  s <- sqrt(diag(th_y_list$Sigma))
  th_y_vect <- c(th_y_vect,s)
  if(cdep) {
    # Then need to add z
    # The leading, intragroup components of z come from groups with two or more
    # members (non-singleton groups)
#    numIntra <- get_num_nonsing_groups(modSpec$cdepGroups)
    nsGroups <-  get_non_singleton_groups(modSpec$cdepGroups) # The non-singleton groups
    zns <- rep(NA,length(nsGroups))
    for(n in 1:length(nsGroups)) {
      g <- nsGroups[n]
      matches <- which(modSpec$cdepGroups == g)
      i1 <- matches[1]
      i2 <- matches[2]
      zns[n] <- th_y_list$Sigma[i1,i2]/s[i1]/s[i2]
    }
    zcross <- rep(NA,choose(numGroups,2))
    for(n in 1:choose(numGroups,2)) {
      e <- elem(n-1,numGroups,2) + 1
      g1 <- e[1]
      g2 <- e[2]
      matches1 <- which(modSpec$cdepGroups == g1)
      matches2 <- which(modSpec$cdepGroups == g2)
      i1 <- matches1[1]
      i2 <- matches2[1]
      zcross[n] <- th_y_list$Sigma[i1,i2]/s[i1]/s[i2]
    }
    th_y_vect <- c(th_y_vect,zns,zcross)
  }

  if(hetero) {
    # Get the first element for each group
    numGroups <- length(unique(modSpec$hetGroups))
    kappa <- rep(NA,numGroups)
    for(g in 1:numGroups) {
      kappa[g] <- th_y_list$kappa[which(th_y_list$modSpec$hetGroups == g)[1]]
    }
    th_y_vect <- c(th_y_vect,kappa)
  }

  return(th_y_vect)
}




#' @export
theta_y_constr2unconstr <- function(th_y_vect,modSpec) {
  # Constraints:
  #
  # rho			positive
  # tau[1]		unconstrained
  # tau[m+1]-tau[m]	positive [m > 1]
  # a			positive 
  # r			positive 
  # b			unconstrained
  # s			positive
  # z			-1 to 1
  # kappa		positive
  #
  # To constrain the variable v to be positive:
  #
  # v    = exp(vbar)
  # vbar = log(v)
  #
  # To constrain the variable v to be -1 to 1:
  #
  # v    = -1 + 2/(1 + exp(-vbar))
  # vbar = logit((z+1)/2)
  check_model(modSpec)

  hetero <- is_hetero(modSpec)
  cdep   <- is_cdep  (modSpec)

  J <- get_J(modSpec)
  K <- get_K(modSpec)


  # rho should be positive 
  if(J > 0) {
    rho_bar <- log(th_y_vect[get_var_index('rho',modSpec)])
  } else {
    rho_bar <- c()
  }

  # tau_1 is unconstrained. Successive differences should be positive
  tau_bar <- list()
  if(J > 0) {
    for(j in 1:modSpec$J) {
      tau_j <- th_y_vect[get_var_index('tau',modSpec,j=j)]
      M_j   <- length(tau_j)
      tau_bar[[j]] <- tau_j[1]
      if(M_j > 1) {
        tau_bar[[j]] <- c(tau_bar[[j]],log(tau_j[2:M_j] - tau_j[1:(M_j-1)]))
      }
    }
    tau_bar <- unlist(tau_bar)
  } else {
    tau_bar <- c()
  }

  # a and r should be positive 
  if(K > 0) {
    a_bar <- log(th_y_vect[get_var_index('a',modSpec)])
    r_bar <- log(th_y_vect[get_var_index('r',modSpec)])
    b <- th_y_vect[get_var_index('b',modSpec)]
  } else {
    a_bar <- c()
    r_bar <- c()
    b <- c()
  }

  # s should be positive 
  s_bar <- log(th_y_vect[get_var_index('s',modSpec)])

  # z should be between -1 and 1
  if(cdep) {
    z_bar <- gtools::logit((th_y_vect[get_var_index('z',modSpec)]+1)/2)
  } else {
    z_bar <- c()
  }
  
  # kappa should be positive
  if(hetero) {
    kappa_bar <- log(th_y_vect[get_var_index('kappa',modSpec)])
  } else {
    kappa_bar <- c()
  }
  return(c(rho_bar,tau_bar,a_bar,r_bar,b,s_bar,z_bar,kappa_bar))
}

#' @export
theta_y_unconstr2constr <- function(th_y_vect,modSpec) {
  # Constraints:
  #
  # rho			positive
  # tau[1]		unconstrained
  # tau[m+1]-tau[m]	positive [m > 1]
  # a			positive 
  # r			positive 
  # b			unconstrained
  # s			positive
  # z			-1 to 1
  # kappa		positive
  #
  # To constrain the variable v to be positive:
  #
  # v    = exp(vbar)
  # vbar = log(v)
  #
  # To constrain the variable v to be -1 to 1:
  #
  # v    = -1 + 2/(1 + exp(-vbar))
  # vbar = logit((z+1)/2)
  check_model(modSpec)

  hetero <- is_hetero(modSpec)
  cdep   <- is_cdep  (modSpec)

  J <- get_J(modSpec)
  K <- get_K(modSpec)


  # rho should be positive 
  if(J > 0) {
    rho <- exp(th_y_vect[get_var_index('rho',modSpec)])
  } else {
    rho <- c()
  }

  # tau_1 is unconstrained. Successive differences should be positive
  tau <- list()
  if(J > 0) {
    for(j in 1:J) {
      tau_bar_j <- th_y_vect[get_var_index('tau',modSpec,j=j)]
      M_j   <- length(tau_bar_j)
      if(M_j == 1) {
        tau[[j]] <- tau_bar_j
      } else {
        tau[[j]] <- tau_bar_j[1] + c(0,cumsum(exp(tau_bar_j[2:M_j])))
      }
    }
    tau <- unlist(tau)
  } else {
    tau <- c()
  }

  # a and r should be positive 
  if(K > 0) {
    a <- exp(th_y_vect[get_var_index('a',modSpec)])
    r <- exp(th_y_vect[get_var_index('r',modSpec)])
    b <- th_y_vect[get_var_index('b',modSpec)]
  } else {
    a <- c()
    r <- c()
    b <- c()
  }

  # s should be positive 
  s <- exp(th_y_vect[get_var_index('s',modSpec)])

  # z should be between -1 and 1
  if(cdep) {
    z <- -1 + 2/(1+exp(-th_y_vect[get_var_index('z',modSpec)]))
  } else {
    z <- c()
  }
  
  # kappa should be positive
  if(hetero) {
    kappa <- exp(th_y_vect[get_var_index('kappa',modSpec)])
  } else {
    kappa <- c()
  }
  return(c(rho,tau,a,r,b,s,z,kappa))
}
