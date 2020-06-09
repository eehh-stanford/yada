#' @title Do a maximum likelihood fit for theta_y
#'
#' @description \code{fit_thy} does a maximum likelihood fit for theta_y
#'
#' @details
#'
#' x is a vector of length n and Y is a matrix of size (J+K) by N, where J is
#' the number of ordinal observations and K the number of continuous
#' observations. modSpec, which stands for model specification, specifies the
#' model. In particular, it is a list with the following fields:
#'
#' J          Number of ordinal observations
#' K          Number of continuous observations
#' M          Number of ordinal categories for each ordinal variable (a vector
#'            of length J)
#' meanSpec   Specification of the parametetric form of the mean response.
#'            Currently, only powLaw, for power law, is supported.
#' hetSpec    Heteroskedasticity specification. Currently, only none (yielding
#'            a homoskedastic model) and linearSd (in which the noise standard
#'            deviation scales linearly with x) are supported. For linearSd,
#'            the variable hetGroups must also be specified.
#' hetGroups  The groups for the heteroskedastic parameters (kappa). A vector
#'            of length J+K. Each variable must be assigned a unique group 
#'            between 1 and G_kappa. Only required if hetSpec is linearSd.
#' cdepSpec   Conditional dependence specification, which equals either indep
#'            (for independent) or dep (for dependent). For dep, the variable
#'            cdepGroups must be specified.
#' cdepGroups The groups for the correlation terms for a conditional dependence
#'            setting of dep. A vector of length J+K. Each variable must be
#'            assigned a unique group between 1 and G_rho. Only required if
#'            cdepSpec is dep.
#'
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param modSpec A list specifying the model to use
#'
# @keywords
#' @export
#'
#' @return The log-likelihood
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
fit_theta_y <- function(x,Y,modSpec,verbose=F) {
  # dep	multi	-	-	-	-	Call bbb
  # indep	single	-	ord	-	Call fitPowLawOrd
  # indep	single	-	cont	-	Call fitPowLaw
  # indep	multi	homo	-	-	recurse
  # indep	multi	hetero	-	Gkap=1	Do a tailored fit
  # indep	multi	hetero	-	Gkap>1	recurse
  check_model(modSpec)
  J <- get_J(modSpec)
  K <- get_K(modSpec)
  cdep   <- is_cdep (modSpec)

  if(cdep) {
    # The model is conditionally depenendent
    stop('Not yet implemented')
  }

  # If this point is reached, the model is conditionally independent
  multi <- is.matrix(Y) # Is this a multi-variable model?
  hetero <- is_hetero(modSpec)

  if(!multi) {
    # Single variable model. Call either fitPowLawOrd or fitPowLaw
    if(J+K != 1) {
      stop('Y is not a matrix yet J and K do not sum to 1')
    }
    if(verbose) {
      if(J == 1) {
        print('Fitting a single-variable, ordinal model')
      } else {
        print('Fitting a single-variable, continuous model')
      }
    }
    # Remove possible missing values
    keep <- !is.na(Y) & !is.na(x)
    if(J == 1) {
      return(fitPowLawOrd(x[keep],Y[keep],hetSpec=modSpec$hetSpec))
    } else {
      return(fitPowLaw(x[keep],Y[keep],hetSpec=modSpec$hetSpec))
    }
  }

  # If this point is reached, the model is conditionally independent and
  # multi-variable.
  if(!hetero) {
    # The model is conditionally independent, multi-variable, and
    # homoskedastic. Subset and recursively call fit_theta_y for each variable.
    if(verbose) {
      print(paste0('Fitting a conditionally independent, multi-variable, homoskedastic model for which J = ',J,' and K = ',K))
    }
    if(J > 0) {
      rho <- rep(NA,J)
      tau <- list()
    } else {
      rho <- c()
      tau <- c()
    }

    if(K > 0) {
      a <- rep(NA,K)
      r <- rep(NA,K)
      b <- rep(NA,K)
    } else {
      a <- c()
      r <- c()
      b <- c()
    }
    s  <- rep(NA,J+K)

    if(J > 0) {
      for(j in 1:J) {
        if(verbose) {
          print(paste('j =',j,'of',J))
        }
        xj <- x
        vj <- Y[j,]
        xj <- xj[!is.na(vj)]
        vj <- vj[!is.na(vj)]
        modSpec_j <- list(meanSpec='powLaw')
        modSpec_j$J <- 1
        modSpec_j$M <- length(unique(vj)) - 1
        modSpec_j$hetSpec <- 'none'
        modSpec_j$cdepSpec <- 'indep'

        th_v <- fit_theta_y(xj,vj,modSpec_j,verbose=verbose)
        th_v_list <- theta_y_vect2list(th_v,modSpec_j)
        rho[j] <- th_v_list$rho
        tau[[j]] <- th_v_list$tau[[1]]
        s[j] <- th_v_list$s
      }
    }

    if(K > 0) {
      for(k in 1:K) {
        if(verbose) {
          print(paste('k =',k,'of',K))
        }
        xk <- x
        vk <- Y[J+k,]
        xk <- xk[!is.na(vk)]
        vk <- vk[!is.na(vk)]
        modSpec_k <- list(meanSpec='powLaw')
        modSpec_k$K <- 1
        modSpec_k$hetSpec <- 'none'
        modSpec_k$cdepSpec <- 'indep'
        th_w <- fit_theta_y(xk,vk,modSpec_k,verbose=verbose)
        th_w_list <- theta_y_vect2list(th_w,modSpec_k)
        a[k] <- th_w_list$a
        r[k] <- th_w_list$r
        b[k] <- th_w_list$b
        s[J+k] <- th_w_list$s
      }      
    }
    return(c(rho,unlist(tau),a,r,b,s))
  }

  # If this point is reached, the model is conditionally independent,
  # multi-variable, and heteroskedastic.

  # Check whether any variables are specified as homoskedastic (NA in hetGroups)
  haveNA <- any(is.na(modSpec$hetGroups))

  if(haveNA) {
    ind_homo <- which(is.na(modSpec$hetGroups))
    J_homo <- sum(ind_homo <= modSpec$J)
    K_homo <- sum(ind_homo > modSpec$J)

    ind_hetero <- which(!is.na(modSpec$hetGroups))
    J_hetero <- sum(ind_hetero <= modSpec$J)
    K_hetero <- sum(ind_hetero > modSpec$J)
    if(verbose) {
      print(paste0('Fitting a conditionally independent, multi-variable, heteroskedastic model for which J = ',J,' and K = ',K))
      print(paste0('Of these variables, ',J_homo,' ordinal variables and ',K_homo,' continuous variables are specially specified as homoskedastic'))
    }
    # Subset and solve separately for the homoskedastic and heteroskedastic variables
    modSpec_homo <- list(meanSpec='powLaw')
    modSpec_homo$cdepSpec <- 'indep'
    modSpec_homo$J <- J_homo
    modSpec_homo$K <- K_homo
    if(J_homo > 0) {
      modSpec_homo$M <- modSpec$M[ind_homo]
    }
    modSpec_homo$hetSpec <- 'none'
    th_y_homo <- fit_theta_y(x,Y[ind_homo,],modSpec_homo,verbose)
    th_y_list_homo <- theta_y_vect2list(th_y_homo,modSpec_homo)

    modSpec_hetero <- list(meanSpec='powLaw')
    modSpec_hetero$cdepSpec <- 'indep'
    modSpec_hetero$J <- J_hetero
    modSpec_hetero$K <- K_hetero
    modSpec_hetero$hetSpec <- modSpec$hetSpec
    modSpec_hetero$hetGroups <- modSpec$hetGroups[ind_hetero]
    if(J_hetero > 0) {
      modSpec_hetero$M <- modSpec$M[ind_hetero]
    }
    modSpec_hetero$hetSpec <- modSpec$hetSpec
    modSpec_hetero$hetGroups <- modSpec$hetGroups[ind_hetero]
    th_y_hetero <- fit_theta_y(x,Y[ind_hetero,],modSpec_hetero,verbose)
    th_y_list_hetero <- theta_y_vect2list(th_y_hetero,modSpec_hetero)
    
    # Combine the results
    th_y_list_full <- list(modSpec=modSpec)
    s <- rep(NA,J+K)

    if(J > 0) {
      rho <- rep(NA,J)
      tau <- list()
      for(j in 1:J) {
        if(is.na(modSpec$hetGroups[j])) {
          # Then j is a homoskedastic variable
          # Index in the homoskedastic fit:
          j_homo <- sum(is.na(modSpec$hetGroups[1:j]))
          rho[j]   <- th_y_list_homo$rho[j_homo]
          tau[[j]] <- th_y_list_homo$tau[[j_homo]]
          s[j]     <- th_y_list_homo$s[j_homo]
        } else {
          # Then j is a heteroskedastic variable
          j_hetero <- sum(!is.na(modSpec$hetGroups[1:j]))
          rho[j]   <- th_y_list_hetero$rho[j_hetero]
          tau[[j]] <- th_y_list_hetero$tau[[j_hetero]]
          s[j]     <- th_y_list_hetero$s[j_hetero]
        }
      }
     th_y_list_full$rho <- rho
     th_y_list_full$tau <- tau
    }

    if(K > 0) {
      a <- rep(NA,K)
      r <- rep(NA,K)
      b <- rep(NA,K)
      for(k in 1:K) {
        if(is.na(modSpec$hetGroups[J+k])) {
          # Then k is a homoskedastic variable
          # Index in the homoskedastic fit:
          k_homo <- sum(is.na(modSpec$hetGroups[1:(J+k)]))
          a[k]   <- th_y_list_homo$a[k_homo]
          r[k]   <- th_y_list_homo$r[k_homo]
          b[k]   <- th_y_list_homo$b[k_homo]
          s[J+k] <- th_y_list_homo$s[k_homo]
        } else {
          # Then j is a heteroskedastic variable
          k_hetero <- sum(!is.na(modSpec$hetGroups[1:(J+k)]))
          a[k]   <- th_y_list_hetero$a[k_hetero]
          r[k]   <- th_y_list_hetero$r[k_hetero]
          b[k]   <- th_y_list_hetero$b[k_hetero]
          s[J+k] <- th_y_list_hetero$s[k_hetero]
        }
      }
      th_y_list_full$a <- a
      th_y_list_full$r <- r
      th_y_list_full$b <- b
    }
    th_y_list_full$s <- s
    th_y_list_full$kappa <- th_y_list_hetero$kappa
    if(modSpec$hetSpec == 'sd_pow') {
      th_y_list_full$lambda <- th_y_list_hetero$lambda
    }
    return(theta_y_list2vect(th_y_list_full))
  }

  # If this point is reached, the model is conditionally independent,
  # multi-variable, heteroskedastic, and contains no special-case homoskedastic
  # variables.
  Gkappa <- get_Gkappa(modSpec)

  if(Gkappa == 1) {
    # The model is conditionally independent, multi-variable, heteroskedastic,
    # and Gkappa = 1. Do a tailored fit.
   if(verbose) {
      print(paste0('Fitting a conditionally independent, multi-variable, heteroskedastic model with a single heteroskedastic parameter'))
   }
    modSpec0 <- list(meanSpec='powLaw')
    modSpec0$cdepSpec <- 'indep'
    modSpec0$hetSpec  <- 'none'
    modSpec0$J        <- modSpec$J
    modSpec0$K        <- modSpec$K
    th_y0 <- fit_theta_y(x,Y,modSpec0,verbose=F)
    th_y0 <- c(th_y0,0.0001) # Add a small amount of heteroskedasticity
    th_y_bar0 <- theta_y_constr2unconstr(th_y0,modSpec)

    # Handle possible missing values
    keep <- colSums(is.na(Y)) < ncol(Y)

    useHjk <- F
    if(useHjk) {
      optimControl <- list(info=verbose)
      fit <- dfoptim::hjk(th_y_bar0,powLawMixNegLogLik,control=optimControl,x=x[keep],Y=Y[,keep],modSpec=modSpec,transformVar=T)
    } else {
      optimControl <- list(reltol=1e-12,maxit=10000000,ndeps=rep(1e-8,length(th_y_bar0)))
      if(verbose) {
        optimControl$trace <- 100
      }
      fit <- optim(th_y_bar0,powLawMixNegLogLik,control=optimControl,x=x[keep],Y=Y[,keep],modSpec=modSpec,transformVar=T,method='BFGS')
    }
  

    th_y <- theta_y_unconstr2constr(fit$par,modSpec)
    return(th_y)
  }

  # If this point is reached, the model is conditionally independent,
  # multi-variable, heteroskedastic, and Gkappa > 1. Subset the problem and
  # recursively call fit_theta_y.

  if(verbose) {
    print(paste0('Fitting a conditionally independent, multi-variable, heteroskedastic model with ',Gkappa,' distinct heteroskedastic parameters'))
  }
  if(J > 0) {
    rho <- rep(NA,J)
    tau <- list()
  } else {
    rho <- c()
    tau <- c()
  }

  if(K > 0) {
    a <- rep(NA,K)
    r <- rep(NA,K)
    b <- rep(NA,K)
  } else {
    a <- list()
    r <- list()
    b <- list()
  }
  s  <- rep(NA,J+K)
  kappa <- rep(NA,Gkappa)
  if(modSpec$hetSpec == 'sd_pow') {
    lambda <- rep(NA,Gkappa)
  } else {
    lambda <- c()
  }

  for(g in 1:Gkappa) {
    if(verbose) {
      print(paste('g =',g,'of',Gkappa))
    }
    # For subsetting variables in this group:
    indg = which(modSpec$hetGroups == g)
    modSpec_g <- list(meanSpec='powLaw')
    J_g <- sum(indg <= J)
    K_g <- sum(indg  > J)
    if(J_g > 0) {
      modSpec_g$J <- J_g
      indg_ord <- indg[indg <= J]
      modSpec_g$M <- modSpec$M[indg[indg <= J]]
    }
    if(K_g > 0) {
      modSpec_g$K <- K_g
      indg_cont <- indg[indg > J]
    }
    modSpec_g$cdepSpec <- 'indep'
    modSpec_g$hetSpec <- modSpec$hetSpec
    modSpec_g$hetGroups <- rep(1,length(indg))
    th_y <- fit_theta_y(x,Y[indg,],modSpec_g,verbose=verbose)

    if(J_g > 0) {
      rho[indg_ord] <- th_y[get_var_index('rho',modSpec_g)]
      tau[[indg_ord]] <- th_y[get_var_index('tau',modSpec_g)]
    }

    if(K_g > 0) {
      a[indg_cont-J] <- th_y[get_var_index('a',modSpec_g)]
      r[indg_cont-J] <- th_y[get_var_index('r',modSpec_g)]
      b[indg_cont-J] <- th_y[get_var_index('b',modSpec_g)]
    }
    s[indg] <- th_y[get_var_index('s',modSpec_g)]
    kappa[g] <- th_y[get_var_index('kappa',modSpec_g)]
    if(modSpec$hetSpec == 'sd_pow') {
      lambda[g] <- th_y[get_var_index('lambda',modSpec_g)]
    }
  }
  return(c(rho,unlist(tau),a,r,b,s,kappa,lambda))
}
