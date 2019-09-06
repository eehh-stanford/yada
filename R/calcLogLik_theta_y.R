#' @title Calculate the log-likelihood for theta_y in the mixed cumulative probit model
#'
#' @description \code{calcLogLik_theta_y} calculates the log-likelihood for theta_y (input as a list) given observed data.
#'
# @details
#'
#' @param theta_y_list theta_y as a list
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param hp Hyperparameters
#'
# @keywords
#' @export
#'
#' @return The log-likelihood
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
calcLogLik_theta_y <- function(theta_y_list,x,Y,hp) {
  '%dopar%' <- foreach::'%dopar%'
  integInfo <- getIntegInfo_theta_y(theta_y_list,Y)
  if(is.matrix(Y)) {
    # For more than one variable, call calcLogLik_theta_y for each observation in parallel
    logLikVect <- foreach::foreach(n=1:ncol(Y), .combine=cbind) %dopar% {
      logLik <- calcLogLik_theta_y(theta_y_list,x[n],Y[,n],hp)
    }
    return(sum(logLikVect))
  } else {
    # The calculation for one observation
    if('gamma' %in% names(theta_y_list)) {
      # using correlations
      covMat <- (1+x)^(1-theta_y_list$gamma) * theta_y_list$Sigma
    } else {
      # not using correlations
      covMat <- theta_y_list$Sigma
    }
    intAll <- all(integInfo$doIntegral)
    # If all variables are integrated, the conditional calculations are not needed
    if(intAll) {
      lo <- integInfo$limArray[,1]
      hi <- integInfo$limArray[,2]
      p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=calc_theta_y_means(x,theta_y_list),sigma=covMat)
      logLik <- log(as.numeric(p))
    } else {
      dep <- which(integInfo$doIntegral) # dep for dependent (conditioned on known variables)
      giv <- which(!integInfo$doIntegral) # giv for given (variables conditioned on)
      lo <- integInfo$limArray[dep,1]
      hi <- integInfo$limArray[dep,2]
      if(length(dep) > 0) {
        # Integral needed
        condNorm <- condMVNorm::condMVN(mean=calc_theta_y_means(x,theta_y_list),sigma=covMat, dependent=dep, given=giv,X.given=Y[giv])
        p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=condNorm$condMean,sigma=condNorm$condVar) # The integral
        logLik <- log(as.numeric(p))
      } else {
        # No integral needed
        logLik <- 0
      }

      if(length(giv) == 1) { # The contribution of the point calculation integral
        logLik <- logLik + log(dnorm(Y[giv],mean=calc_theta_y_means(x,theta_y_list,giv),sd=sqrt(covMat[giv,giv])))
      } else {
        logLik <- logLik + mvtnorm::dmvnorm(Y[giv],mean=calc_theta_y_means(x,theta_y_list,giv),sigma=covMat[giv,giv],log=T)
      }
    }
    return(logLik)
  }
}

# A wrapper function to calculate the mean for each variable
calc_theta_y_means <- function(xScalar,theta_y_list,giv=NA) {
  haveOrd <- 'rho' %in% names(theta_y_list)
  haveCont <- 'r' %in% names(theta_y_list)
  validModels <- tolower(c('GenCRRA_corr_heterosk','GenCRRA_corr_homosk','GenCRRA_uncorr_heterosk','GenCRRA_uncorr_homosk'))

  if(tolower(theta_y_list$paramModel) %in% validModels) {
    if(haveOrd) {
      rho <- theta_y_list$rho
    }

    if(haveCont) {
      a <- theta_y_list$a
      r <- theta_y_list$r
      b <- theta_y_list$b
    }

    if(!any(is.na(giv))) {
      if(haveOrd && haveCont) {
        J <- length(rho)
        givOrd  <- giv[giv <= J]
        givCont <- giv[giv > J]
        rho <- rho[givOrd]
        a <- a[givCont-J]
        r <- r[givCont-J]
        b <- b[givCont-J]
      } else if(haveOrd && !haveCont) {
        givOrd  <- giv
        rho <- rho[givOrd]
      } else {
        # !haveOrd && haveCont
        givCont  <- giv
        a <- a[givCont]
        r <- r[givCont]
        b <- b[givCont]
      }
    }
    if(haveOrd && haveCont) {
      meanVectOrd  <- xScalar^(1-rho)
      meanVectCont <- a*xScalar^(1-r) + b
      return(c(meanVectOrd,meanVectCont))
    } else if(haveOrd && !haveCont) {
      meanVectOrd  <- xScalar^(1-rho)
      return(meanVectOrd)
    } else {
      # !haveOrd && haveCont
      meanVectCont <- a*xScalar^(1-r) + b
      return(meanVectCont)
    }

  } else {
    stop(paste('Unsupported parametric model specification:',hp$paramModel))
  }
}
