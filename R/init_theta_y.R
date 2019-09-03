#' @title Initialize parameter vector for the mixed cumulative probit model
#'
#' @description \code{init_theta_y} initializes the parameter vector, theta_y (e.g., for maximum likelihood estimation or adaptive Metropolis sampling). Individual variable fits are done by calling fitGenCrraOrd and fitGenCrra and all correlations are set to zero (off diagonal terms in Sigma).
#'
# @details
#'
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param J Number of latent variables
#' @param corr Whether to include correlations across latent variables [default True]
#' @param heter0 Whether to allow heteroskedasticity [default True]
#'
# @examples
#'
#' @return theta_y as a list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
#' @export

init_theta_y <- function(x,Y,J,corr=T,hetero=T) {
  K <- nrow(Y) - J # Remaining variables are continuous
  
  if(J > 0) {
    rho <- rep(NA,J)
    sig_ord <- rep(NA,J)
    tau <- list()
    for(j in 1:J) {
      param <- fitGenCrraOrd(x,Y[j,])
      rho[j] <- param[1]
      sig_ord[j] <- param[2]
      tau[[j]] <- param[3:length(param)]
    }
  }

  if(K > 0) {
    a <- rep(NA,K)
    r <- rep(NA,K)
    b <- rep(NA,K)
    sig_cont <- rep(NA,K)
    for(k in 1:K) {
      param <- fitGenCrra(x,Y[J+k,],fitSig=T)
      a[k] <- param[1]
      r[k] <- param[2]
      b[k] <- param[3]
      sig_cont[k] <- param[4]
    }
  }

  if(corr) {
    if(hetero) {
      paramModel <- 'GenCRRA_corr_heterosk' # correlations and heteroskedasticity
    } else {
      paramModel <- 'GenCRRA_corr_homosk' # correlations and homoskedasticity
    }
  } else {
    if(hetero) {
      paramModel <- 'GenCRRA_uncorr_heterosk' # no correlations and heteroskedasticity
    } else {
      paramModel <- 'GenCRRA_uncorr_homosk' # no correlations and homoskedasticity
    }
  }

  theta_y_list <- list(paramModel=paramModel,rho=rho,a=a,r=r,b=b,tau=tau)
  if(J > 0 && K > 0) {
    theta_y_list$Sigma <- diag(c(sig_ord,sig_cont)^2)
  } else if(J > 0 && K == 0) {
    theta_y_list$Sigma <- diag(sig_ord^2)
  } else {
    # J == 0 && K > 0
    theta_y_list$Sigma <- diag(sig_cont^2)
  }

  if(hetero) {
    theta_y_list$gamma <- 0.5
  }
  return(theta_y_list)
}
