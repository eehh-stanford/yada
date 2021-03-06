#' @title Sample theta_y in the mixed cumulative probit model
#'
#' @description \code{sample_theta_y} uses the adaptive Metropolis algorithm of Haario et al. 2001, https://projecteuclid.org/euclid.bj/1080222083, to sample the vector theta_y for the cumulative probit model.
#'
# @details
#'
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param hp Hyperparameters
#' @param numSamp Number of samples
#' @param burnIn Number of burn-in samples
#' @param verbose Whether to print out information as the run proceeds [optional]
#' @param start Starting value of theta_y_list for sampling [optional]
#' @param fileName File name for saving [optional]
#' @param savePeriod Period (of samples) for saving [optiona; default 1000]
#' @param known Known values, a list with the field theta_y, which is also a list [optional]
#' @param varNames Variable names for saving to file [optional]
#'
#' @return theta_y as a list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
# @references

#' @export
sample_theta_y <- function(x,Y,hp,numSamp,burnIn,verbose=T,start=NA,fileName=NA,savePeriod=1000,known=NA,varNames=NA) {
  # Initialize theta_y and get ready for sampling
  if(any(is.na(start))) {
    theta_y_t_list <- init_theta_y(x,Y,hp$J)
  } else {
    theta_y_t_list <- start
  }

  theta_y_t <- theta_yList2Vect(theta_y_t_list)
  logLik_t <- calcLogLik_theta_y(theta_y_t_list,x,Y,hp)

  theta_y_0 <- theta_y_t

  s_d <- (2.4)^2 / length(theta_y_t)
  epsilon <- 1e-12
  I_d <- diag(length(theta_y_t))
  mean_theta_y <- theta_y_t

  c0 <- 1e-6
  C0 <- diag(c(rep(c0,length(theta_y_t))))
  thetaList <- list()
  logLikVect <- vector()
  for(tt in 1:(burnIn+numSamp)) {
    if(tt <= burnIn) {
      if(verbose) {
        print('---- [burn in]')
      }
      theta_y_tp1 <- MASS::mvrnorm(1,mu=theta_y_t,Sigma=C0)
    } else {
      if(verbose) {
        print('----')
      }
      C_t <- s_d * (cov_theta_y + epsilon*I_d)
      theta_y_tp1 <- MASS::mvrnorm(1,mu=theta_y_t,Sigma=C_t)
    }
    if(verbose) {
      print(tt)
    }
    theta_y_tp1_list <- theta_yVect2List(theta_y_tp1,hp)
    
    if(!theta_yIsValid(theta_y_tp1_list,hp)) {
      # Reject the new sample
      #theta_tp1 <- theta_t
      accept <- F 
      print('**')
      print('**')
      print('**')
      print('**')
    } else {
      # Calculate the acceptance ratio
      # By construction, the proposal distribution is symmetric. Hence:
      # A = p(theta_tp1|alpha,D) / p(theta_t|alpha,D)
      # A = p(D|theta_tp1) * p(theta_tp1|alpha) / 
      #     p(D|theta_t) / p(theta_t|alpha) 
      logLik_tp1 <- calcLogLik_theta_y(theta_y_tp1_list,x,Y,hp)
      #logPrior_kp1 <- calcLogPrior(theta_kp1_list,hp)
    
      #logA <- logPrior_kp1 + logLik_kp1 - logPrior_k - logLik_k
      logA <- logLik_tp1 - logLik_t
      if(verbose) {
        if(tt <= burnIn) {
          print(C0[1,1])
          print(C0[100,100])
        } else {
          print(C_t[1,1])
          print(C_t[100,100])
        }
        print(logA)
        print(max(logLik_t,logLik_tp1))
      }
      if(is.finite(logA)) {
        logA <- min(0,logA)
        accept <- log(runif(1)) < logA
      } else {
        accept <- F
      }
    }
    if(!accept) {
      theta_y_tp1 <- theta_y_t
      theta_y_tp1_list <- theta_y_t_list
      #logPrior_kp1 <- logPrior_k
      logLik_tp1 <- logLik_t
    }
    thetaList[[tt]] <- theta_y_tp1_list
    logLikVect[tt] <- logLik_tp1


    # Get ready for next sample
    theta_y_t <- theta_y_tp1
    theta_y_t_list <- theta_y_tp1_list
    #logPrior_k <- logPrior_kp1
    logLik_t <- logLik_tp1
    mean_theta_y <- mean_theta_y * tt/(tt+1) + theta_y_t/(tt+1)
    if(tt == 1) {
      B_t <- tcrossprod(theta_y_0,theta_y_0) + tcrossprod(theta_y_t,theta_y_t)
    } else {
      B_t <- B_t*(tt-1)/tt + tcrossprod(theta_y_t,theta_y_t)/tt
    }
    cov_theta_y <- B_t - (tt+1)*tcrossprod(mean_theta_y,mean_theta_y)/tt
    if(verbose) {
      print(accept)
    }
    if(!is.na(fileName)) {
      if(tt %% savePeriod == 0 ) {
        save_theta_y_sample(fileName,thetaList,logLikVect,x,Y,hp,varNames,known)
      }
    }
  }
  return(list(theta_yList=thetaList,logLikVect=logLikVect))
}
