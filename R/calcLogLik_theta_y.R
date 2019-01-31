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
    logLikVect <- rep(0,length(x))
    logLikVect <- foreach::foreach(n=1:ncol(Y), .combine=cbind) %dopar% {
      logLik <- calcLogLik_theta_y(theta_y_list,x[n],Y[,n],hp)
    }
    return(sum(logLikVect))
  } else {
    # The calculation for one observation
    intAll <- all(integInfo$doIntegral)
    # If all variables are integrated, the conditional calculations are not needed
    if(intAll) {
      lo <- integInfo$limArray[,1]
      hi <- integInfo$limArray[,2]
      p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=calc_theta_y_means(x,theta_y_list),sigma=theta_y_list$Sigma)
      logLik <- log(as.numeric(p))
    } else {
      dep <- which(integInfo$doIntegral) # dep for dependent (conditioned on known variables)
      giv <- which(!integInfo$doIntegral) # giv for given (variables conditioned on)
      lo <- integInfo$limArray[dep,1]
      hi <- integInfo$limArray[dep,2]
      condNorm <- condMVNorm::condMVN(mean=calc_theta_y_means(x,theta_y_list),sigma=theta_y_list$Sigma, dependent=dep, given=giv,X.given=Y[giv])
      p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=condNorm$condMean,sigma=condNorm$condVar) # The integral
      logLik <- log(as.numeric(p))
      if(length(giv) == 1) { # The contribution if no integral
        logLik <- logLik + log(dnorm(Y[giv],mean=calc_theta_y_means(x,theta_y_list,giv),sd=sqrt(theta_y_list$Sigma[giv,giv])))
      } else {
        logLik <- logLik + mvtnorm::dmvnorm(Y[giv],mean=calc_theta_y_means(x,theta_y_list,giv),sigma=theta_y_list$Sigma[giv,giv],log=T)
      }
    }
  }
}

# A wrapper function to calculate the mean for each variable
calc_theta_y_means <- function(xScalar,theta_y_list,giv=NA) {
  if(tolower(theta_y_list$paramModel) == 'gencrra') {
    alpha <- theta_y_list$alpha
    rho <- theta_y_list$rho
    a <- theta_y_list$a
    r <- theta_y_list$r
    b <- theta_y_list$b
    if(!any(is.na(giv))) {
      J <- length(alpha)
      givOrd  <- giv[giv <= J]
      givCont <- giv[giv > J]
      alpha <- alpha[givOrd]
      rho <- rho[givOrd]
      a <- a[givCont-J]
      r <- r[givCont-J]
      b <- b[givCont-J]
    }

    meanVectOrd  <- alpha*xScalar^(1-rho)
    meanVectCont <- a*xScalar^(1-r) + b
    return(c(meanVectOrd,meanVectCont))
  } else {
    stop(paste('Unsupported parametric model specification:',hp$paramModel))
  }
}
