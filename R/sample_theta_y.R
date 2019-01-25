#' @title Sample theta_y in the mixed cumulative probit model
#'
#' @description \code{sample_theta_y} uses the adaptive Metropolis algorithm of Haario et al. 2001, https://projecteuclid.org/euclid.bj/1080222083, to sample the vector theta_y for the cumulative probit model.
#'
# @details
#'
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param J Number of latent variables
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return theta_y as a list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
# @references


sample_theta_y <- function(x,Y,hp,numSamp,burnIn,verbose=T) {
  # Initialize theta_y and get ready for sampling
  theta_y_t_list <- init_theta_y(x,Y,hp$J)
  theta_y_t <- theta_yList2Vect(theta_y_t_list)
  logLik_t <- calcLogLik_theta_y(theta_y_t_list,x,Y,hp)

  theta_y_0 <- theta_y_t

  s_d <- (2.4)^2 / length(theta_y_t)
  epsilon <- 1e-10
  I_d <- diag(length(theta_y_t))
  mean_theta_y <- theta_y_t

  c0 <- 1e-6
  C0 <- diag(c(rep(c0,length(theta_y_t))))
  thetaList <- list()
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
#    if(tt %% summaryPeriod == 1 ) {
#      known <- list(theta=theta,Ystar=Ystar)
#      saveMH('mh_save.rds',thetaList,Y,x,hp,vars,known)
#      mh <- readRDS('mh_save.rds')
#      summarizeResults(mh,floor(tt/1000))
#    }
  }

  return(thetaList)
}
