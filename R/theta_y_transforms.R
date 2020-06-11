#' @title Transformation functions for theta_y
#'
#' @description A number of transformation functions for theta_y -- for example, converting from a list to vector representation, and vice versa.
#'
#' @details For for variable definitions see powLawMix.R. For every ordinal variable, rho, tau, and s are uniquely specified. For every continuous variable, a, r, b, and s are uniquely specified. For all variables, kappa (the heteroskedastic parameter) and z (the correlation parameter) can be grouped or set to zero using the vectors heteroSpec and cdepSpec in the model specification, modSpec.
#'
#' For example, consider a model with two ordinal and two continuous variables with the following model specification (this can be obtained by calling the function yada::get_example_modSpec):
#'
#' modSpec <- list(meanSpec = 'powLaw')
#' modSpec$J <- 2
#' modSpec$K <- 2
#' modSpec$M <- c(1,2)
#' modSpec$hetSpec <- 'sd_x'
#' modSpec$hetGroups <- c(1,NA,2,1)
#' modSpec$cdepSpec <- 'dep'
#' modSpec$cdepGroups <- c(1,1,2,2)
#'
#' To obtain and check this model use:
#'
#' modSpec <- get_example_modSpec()
#' check_model(modSpec)
#'
#' For this model specification, the full heteroskedastic parameter vector is:
#'
#' kappa_full = [kappa1; 0; kappa2; kappa1]
#'
#' The full covariance matrix is:
#'
#' Sigma_full = [   s1*s1, z1*s1*s2, z3*s1*s3, z3*s1*s4;
#'               z1*s2*s1,    s2*s2, z3*s2*s3, z3*s2*s4;
#'               z3*s3*s1, z3*s3*s2,    s3*s3, z2*s3*s4;
#'               z3*s4*s1, z3*s4*s2, z2*s4*s3,    s4*s4]
#'
#' The specification of kappa is relatively straightforward. For variable i, kappa_full[i] = 0 if hetSpec[i] is NA. Otherwise, kappa_full[i] is kappa[hetSpec[i]].
#'
#' The specification of the covariance matrix correlation terms is more involved. If cdpeGroups[i] is NA then all correlation terms associated with that variable are 0. Beyond this, there are two components to the reduced correlation vector z: intra-group correlations and inter-group correlations.
#'
#' The intra-group correlations are for non-singleton groups -- that is, groups with more than one entry in cdepGroups. Let Gns be the number of such non-singleton groups. The first entries of the reduced z are these correlations. In the example above, they are zns = [z1; z2].
#'
#' The inter-group correlations link all groups, whether singleton or not. Let Gz be the number of unique groups of all types. In the example above, Gz = 2. The number of unique inter-group correlations is choose(Gz,2). In the example above, there is one cross-term, zcross = [z3], and z = [zns; zcros] = [z1; z2; z3].
#'
#' The full specification for z (a vector with choose(J+K,2) terms) is
#'
#' zfull = [z1; z3; z3; z3; z3; z2]
#'
#' These correspond to the correlation terms in the upper-right part of Sigma_full, "unwrapped" by rows. Identically, the components of zfull can be thought of as the lexically ordered unique variable pairs {11,12,13,14,23,24,34}; in each pair, the first element gives the row index in Sigma_full and the second element the column index.
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
get_example_modSpec <- function() {
  modSpec <- list(meanSpec = 'powLaw')
  modSpec$J <- 2
  modSpec$K <- 2
  modSpec$M <- c(1,2)
  modSpec$hetSpec <- 'sd_x'
  modSpec$hetGroups <- c(1,NA,2,1)
  modSpec$cdepSpec <- 'dep'
  modSpec$cdepGroups <- c(1,1,2,2)
  return(modSpec)
}

#' @export
check_model <- function(modSpec) {

  if(tolower(modSpec$meanSpec) != 'powlaw') {
    stop(paste('Unrecognized specification for the mean,',modSpec$meanSpec))
  }

  J <- get_J(modSpec)
  K <- get_K(modSpec)

  # Check M
  if(J > 0) {
    if(!('M' %in% names(modSpec))) {
      stop('M should be a field in modSpec if J > 0')
    }

    if(any(is.na(modSpec$M))) {
      stop('No element of M should be NA')
    }

    if(J != length(modSpec$M)) {
      stop(paste0('J [',J,'] should equal length of M [',length(modSpec$M),']'))
    }
  }

  if(J + K == 1) {
    # Single variable models should either have no specification of the
    # conditional dependence or have it be independent
    if('cdepSpec' %in% names(modSpec)) {
      if(tolower(modSpec$cdepSpec) == 'dep') {
        stop('Model is single-variable, but conditional dependence is specified')
      }
    }
  } else {
    # Multiple variable models must have a specification of the conditional
    # dependence, and it must be 'indep' or 'dep'
    if( !('cdepSpec' %in% names(modSpec))) {
      stop('Model is multi-variable, yet cdepSpec is not in modSpec')
    } else {
      if( !(modSpec$cdepSpec %in% c('indep','dep')) ) {
        stop(paste0('cdepSpec should be indep or dep, not ',modSpec$cdepSpec))
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

  # hetSpec must always be specified
  if( !('hetSpec' %in% names(modSpec)) ) {
    stop('hetSpec must be specified')
  }

  # Reject models with an unrecognized hetSpec
  if( !(modSpec$hetSpec %in% c('none','sd_x','sd_resp')) ) {
    stop(paste0('hetSpec = ',modSpec$hetSpec,' is an unrecognized specification'))
  }

  # Reject models specified as heteroskedastic for which hetGroups is
  # mis-specified.
  if(modSpec$hetSpec != 'none') {
    # hetGroups must be specified
    if( !('hetGroups' %in% names(modSpec)) ) {
      stop('Model is heteroskedastic, but hetGroups not given')
    }

    # hetGroups must must be the same length as the number of variables
    N <- length(modSpec$hetGroups)
    if(N != J + K) {
      stop(paste('Length of hetGroups =',N,'but should be J+K=',J+K))
    }

    # hetGroups should be unique, sequential integers starting from 1
    Gkappa <- max(modSpec$hetGroups,na.rm=T)
    uniqueVal <- sort(unique(modSpec$hetGroups[!is.na(modSpec$hetGroups)]))
    if(!all(uniqueVal == 1:Gkappa)) {
      stop('hetGroups should be unique, sequential integers starting from 1')
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
  } else if(tolower(modSpec$hetSpec) == 'sd_x') {
    return(!all(is.na(modSpec$hetGroups)))
  } else if(tolower(modSpec$hetSpec) == 'sd_resp') {
    return(!all(is.na(modSpec$hetGroups)))
  } else {
    stop(paste('Unsupported hetSpec,',modSpec$hetSpec))
  }
}

#' @export
is_cdep <- function(modSpec) {
  # For the special case of a one-variable model, False is returned. Aside from
  # this, the model is conditionally indepenendent if (a) modSpec$cdepSpec is
  # 'indep' or (b) modSpec$cdepSpec is 'dep', but cdepGroups is all NA.

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
get_non_singleton_groups <- function(groupVect) {
  # groupVect is a vector assigning each of its elements to a unique group.
  # Return the non-singleton groups in groupVect -- that is, the groups with
  # more than one member
  #
  # For example:
  # groupVect <- c(1,2,NA,4,1,3,3)
  # print(get_non_singleton_groups(groupVect))
  # [1] 1 3
  counts <- table(groupVect)
  return(sort(as.numeric(names(counts))[counts > 1]))
}

#' @export
get_num_var <- function(varName,modSpec,preceding=F) {
  # Get the number of variables give the model specification for the following
  # variables:
  #
  # rho tau a r b s z kappa
  #
  # If the variable is not part of the model, 0 is returned (rather than
  # throwing an error). This happens, for example, with rho if J=0 or with
  # kappa if the model is homoskedastic.
  #
  # If preceding is TRUE, the number of variables of that precede this variable
  # in the vector is returned.
  check_model(modSpec)

  variables <- c('rho','tau','a','r','b','s','z','kappa')

  if(!(varName %in% variables)) {
    stop(paste('Unrecognized variable name,',varName))
  }

  if(preceding) {
    if(varName == 'rho') {
      return(0)
#    } else if(varName == 'tau') {
#      return(get_num_var('rho',modSpec,preceding=F))
    } else {
      ind <- which(variables == varName) # index in variables
      return(get_num_var(variables[ind-1],modSpec,preceding=T) + get_num_var(variables[ind-1],modSpec,preceding=F))
    }
  }

  if(varName == 'rho') {
    return(get_J(modSpec))
  } else if(varName == 'tau') {
    J <- get_J(modSpec)
    if(J==0) {
      return(0)
    } else {
      return(sum(modSpec$M))
    }
  } else if(varName %in% c('a','r','b')) {
    return(get_K(modSpec))
  } else if(varName == 's') {
    return(get_J(modSpec) + get_K(modSpec))
  } else if(varName == 'z') {
    return(get_z_length(modSpec))
  } else if(varName == 'kappa') {
    return(get_Gkappa(modSpec))
  }
  stop('This point should never be reached')
}

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

  # Handle z specially before handling other variables
  if(varName == 'z') {
    if(!cdep) {
      stop('z requested but model is not conditionally dependent')
    }

    offset <- get_num_var('z',modSpec,preceding=T)
    if(is.na(i1) && is.na(i2)) {
      # Then return all the indices
      return(offset + 1:get_num_var('z',modSpec))
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
    numGroups <- length(unique(modSpec$cdepGroups)[!is.na(unique(modSpec$cdepGroups))])
    if(g1 < g2) {
      index <- elemToIndex(c(g1-1,g2-1),numGroups) + 1
    } else {
      index <- elemToIndex(c(g2-1,g1-1),numGroups) + 1
    }
    return(offset+index)
  }

  # Also handle kappa specially before handling other variables
  if(varName == 'kappa') {
    if(!hetero) {
      stop('kappa requested but model is not heteroskedastic')
    }

    # If modSpec$M is not in M, sum(modSpec$M) evaluates to 0
    offset <- get_num_var('kappa',modSpec,preceding=T)

    if(is.na(j) && is.na(k)) {
      # Then return all the indices
      return(offset + (1:get_num_var('kappa',modSpec)))
    }

    if(!is.na(j)) {
      value <- modSpec$hetGroups[j]
    } else {
      value <- modSpec$hetGroups[J+k]
    }

    # If this variable is not heteroskedastic, return NA
    if(is.na(value)) {
      return(NA)
    } else {
      return(offset+value)
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

  # If neither j nor k is specified, return all indices for the variable
  if(is.na(j) && is.na(k)) {
    return(get_num_var(varName,modSpec,preceding=T) + (1:get_num_var(varName,modSpec)))
  }

  # th_y has ordering th_y = [rho,tau,a,r,b,kappa,s,z]
  offset <- get_num_var(varName,modSpec,preceding=T)
  if(varName == 'rho') {
    return(offset+j)
  } else if(varName == 'tau') {
    if(j > 1) {
      offset <- offset + sum(modSpec$M[1:(j-1)])
    }
    return((offset + 1):(offset+modSpec$M[j]))
  } else if(varName %in% c('a','r','b')) {
    return(offset + k)
  } else if(varName == 's') {
    #         rho   M                a/r/b
    offset <- 1*J + sum(modSpec$M) + 3*K
    if(!is.na(j)) {
      return(offset + j)
    } else {
      return(offset + J + k)
    }
  } else {
    stop(paste('Unrecognized variable',varName))
  }
}

#' @export
theta_y_list2vect <- function(th_y_list) {
  modSpec <- th_y_list$modSpec # Extract to improve code readability
  check_model(modSpec)
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

  th_y_vect <- c(th_y_vect,th_y_list$s)
  if(cdep) {
    # Then need to add z
    th_y_vect <- c(th_y_vect,th_y_list$z)
  }

  if(hetero) {
    # Then need to add kappa
    th_y_vect <- c(th_y_vect,th_y_list$kappa)
  }

  return(th_y_vect)
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

  th_y_list$s <- s

  # For a conditionally dependent model, add z
  if(cdep) {
    th_y_list$z <- th_y_vect[get_var_index('z',modSpec)]
  }

  # For a heteroskedastic model, add kappa
  if(hetero) {
    th_y_list$kappa <- th_y_vect[get_var_index('kappa',modSpec)]
  }

  return(th_y_list)
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

#' @export
get_kappa_full <- function(th_y,modSpec=NA,asMatrix=F) {
  # If modSpec is NA (not input), then the input is th_y_list. If it is given,
  # it is th_y_vect.
  if(!all(is.na(modSpec))) {
    th_y_vect <- th_y
    th_y_list <- theta_y_vect2list(th_y,modSpec)
  } else {
    th_y_list <- th_y
    th_y_vect <- theta_y_list2vect(th_y_list)
  }
  
  check_model(th_y_list$modSpec)

  modSpec <- th_y_list$modSpec
  if(!is_hetero(modSpec)) {
    stop('get_kappa_full called, but model is not heteroskedastic')
  }


  J <- get_J(modSpec)
  K <- get_K(modSpec)
  kappa <- rep(NA,J+K)

  if(J > 0) {
    for(j in 1:J) {
      ind <- get_var_index('kappa',modSpec,j=j)
      if(is.na(ind)) {
        kappa[j] <- 0
      } else {
        kappa[j] <- th_y_vect[ind]
      }
    }
  }

  if(K > 0) {
    for(k in 1:K) {
      ind <- get_var_index('kappa',modSpec,k=k)
      if(is.na(ind)) {
        kappa[J+k] <- 0
      } else {
        kappa[J+k] <- th_y_vect[ind]
      }
    }
  }

  if(!asMatrix) {
    return(kappa)
  } else {
    kappa <- matrix(kappa) # a column vector
    return(kappa %*% t(kappa))
  }
}

#' @export
get_z_full <- function(th_y,modSpec=NA,asMatrix=F) {
  # If modSpec is NA (not input), then the input is th_y_list. If it is given,
  # it is th_y_vect.
  if(!all(is.na(modSpec))) {
    th_y_vect <- th_y
    th_y_list <- theta_y_vect2list(th_y,modSpec)
  } else {
    th_y_list <- th_y
    th_y_vect <- theta_y_list2vect(th_y_list)
  }
  
  check_model(th_y_list$modSpec)

  modSpec <- th_y_list$modSpec
  if(!is_cdep(modSpec)) {
    stop('get_z_full called, but model is not conditionally dependent')
  }


  J <- get_J(modSpec)
  K <- get_K(modSpec)
  z <- rep(NA,choose(J+K,2))

  for(n in 1:choose(J+K,2)) {
    e <- elem(n-1,J+K,2)
    i1 <- e[1] + 1
    i2 <- e[2] + 1
    ind <- get_var_index('z',modSpec,i1=i1,i2=i2)
    if(is.na(ind)) {
      z[n] <- 0
    } else {
      z[n] <- th_y_vect[ind]
    }
  }
  return(z)
}

#' @export
get_Sigma0 <- function(th_y,modSpec=NA) {
  # If modSpec is NA (not input), then the input is th_y_list. If it is given,
  # it is th_y_vect.
  if(!all(is.na(modSpec))) {
    th_y_vect <- th_y
    th_y_list <- theta_y_vect2list(th_y,modSpec)
  } else {
    th_y_list <- th_y
    th_y_vect <- theta_y_list2vect(th_y_list)
  }
  
  check_model(th_y_list$modSpec)

  modSpec <- th_y_list$modSpec
  s <- th_y_list$s
  if(length(s) == 1) {
    return(matrix(s^2,1,1))
  }
  Sigma <- diag(s^2)

  if(!is_cdep(modSpec)) {
    return(Sigma)
  }

  # If this point is reached, the model is conditionally dependent. Add
  # off-diagonal terms.
  J <- get_J(modSpec)
  K <- get_K(modSpec)

  for(i1 in 1:(J+K-1)) {
    g1 <- modSpec$cdepGroups[i1]
    for(i2 in (i1+1):(J+K)) {
      g2 <- modSpec$cdepGroups[i2]
      if(is.na(g1) || is.na(g2)) {
        z12 <- 0
      } else {
        z12 <- th_y_vect[get_var_index('z',modSpec,i1=i1,i2=i2)]
      }
      Sigma[i1,i2] <- s[i1]*s[i2]*z12
      Sigma[i2,i1] <- Sigma[i1,i2]
    }
  }
  return(Sigma)
}

#' @export
get_Sigma <- function(th_y,x,modSpec=NA,transformVar=F) {
  # If modSpec is NA (not input), then the input is th_y_list. If it is given,
  # it is th_y_vect.
  if(all(is.na(modSpec))) {
    modSpec <- th_y$modSpec
    th_y_vect <- theta_y_list2vect(th_y)
  } else {
    th_y_vect <- th_y
  }

  if(transformVar) {
    th_y_vect <- theta_y_unconstr2constr(th_y_vect,modSpec)
  }

  th_y_list <- theta_y_vect2list(th_y_vect,modSpec)
 

  check_model(th_y_list$modSpec)
  J <- get_J(th_y_list$modSpec)
  K <- get_K(th_y_list$modSpec)

  N <- length(x)
  if(N > 1) {
    Sigma <- array(NA,c(N,J+K,J+K))
    for(n in 1:N) {
      Sigma[n,,] <- get_Sigma(th_y_list,x[n])
    }
    return(Sigma)
  }

  # If this point is reached, x is a scalar
  Sigma0 <- get_Sigma0(th_y_list)

  if(!is_hetero(modSpec)) {
    return(Sigma0)
  }
 
  # If this point is reached, x is a scalar and the model is
  # heteroskedastic
  kappa_full <- get_kappa_full(th_y_list)
  heteroTerm <- rep(NA,length(kappa_full)) # length J+K

  if(modSpec$hetSpec == 'sd_x') {
    heteroTerm <- 1 + kappa_full*x
  } else if(modSpec$hetSpec == 'sd_resp') {
    respVect <- c()
    if(J > 0) {
      respVect <- c(respVect,x^th_y_list$rho)
    }
    if(K > 0) {
      respVect <- c(respVect,th_y_list$a*x^th_y_list$r)
    }
    heteroTerm <- 1 + kappa_full*respVect
  } else {
    stop(paste('Unsupported hetSpec,',th_y_list$modSpec$hetSpec))
  }
  # A column vector:
  heteroTerm <- matrix(heteroTerm)
  return((heteroTerm %*% t(heteroTerm) ) * Sigma0)
}
