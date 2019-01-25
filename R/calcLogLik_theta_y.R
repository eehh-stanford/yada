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
  integInfo <- getIntegInfo_theta_y(theta_y_list,Y)

  logLikVect <- rep(0,length(x))

  # Calculate the log-likelihood using a parallel for loop
  logLikVect <- foreach(n=1:ncol(Y), .combine=cbind) %dopar% {
    logLik <- 0
    intAll <- all(integInfo$doIntegral[,n])
    # If all variables are integrated, the conditional calculations are not needed
    if(intAll) {
      lo <- integInfo$limArray[,n,1]
      hi <- integInfo$limArray[,n,2]
      p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=calc_theta_y_means(x[n],theta_y_list),sigma=theta_y_list$Sigma)
      logLik <- logLik + log(as.numeric(p))
    } else {
      dep <- which(integInfo$doIntegral[,n]) # dep for dependent (conditioned on known variables)
      giv <- which(!integInfo$doIntegral[,n]) # giv for given (variables conditioned on)
      lo <- integInfo$limArray[dep,n,1]
      hi <- integInfo$limArray[dep,n,2]
      condNorm <- condMVNorm::condMVN(mean=calc_theta_y_means(x[n],theta_y_list),sigma=theta_y_list$Sigma, dependent=dep, given=giv,X.given=Y[giv,n])
      p <- mvtnorm::pmvnorm(lower=lo, upper=hi,mean=condNorm$condMean,sigma=condNorm$condVar) # The integral
      logLik <- logLik + log(as.numeric(p))
      if(length(giv) == 1) { # The contribution if no integral
        logLik <- logLik + log(dnorm(Y[giv,n],mean=calc_theta_y_means(x[n],theta_y_list,giv),sd=sqrt(theta_y_list$Sigma[giv,giv])))
      } else {
        logLik <- logLik + mvtnorm::dmvnorm(Y[giv,n],mean=calc_theta_y_means(x[n],theta_y_list,giv),sigma=theta_y_list$Sigma[giv,giv],log=T)
      }
    }
    logLik <- logLik
  }
 # print(logLikVect)
  return(sum(logLikVect))
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
