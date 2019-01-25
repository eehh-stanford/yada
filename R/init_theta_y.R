#' @title Initialize parameter vector for the mixed cumulative probit model
#'
#' @description \code{init_theta_y} initializes the parameter vector, theta_y (e.g., for maximum likelihood estimation or adaptive Metropolis sampling). Individual variable fits are done by calling fitGenCrraOrd and fitGenCrra and all correlations are set to zero (off diagonal terms in Sigma).
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


init_theta_y <- function(x,Y,J) {
  K <- nrow(Y) - J # Remaining variables are continuous
  
  alpha <- rep(NA,J)
  rho <- rep(NA,J)
  tau <- list()
  for(j in 1:J) {
    param <- fitGenCrraOrd(x,Y[j,])
    alpha[j] <- param[1]
    rho[j] <- param[2]
    tau[[j]] <- param[3:length(param)]
  }

  a <- rep(NA,hp$K)
  r <- rep(NA,hp$K)
  b <- rep(NA,hp$K)
  for(k in 1:hp$K) {
    param <- fitGenCrra(x,Y[J+k,])
    a[k] <- param[1]
    r[k] <- param[2]
    b[k] <- param[3]
  }

  Sigma <- diag(J+K)
  
  return(list(paramModel='GenCRRA',alpha=alpha,rho=rho,a=a,r=r,b=b,tau=tau,Sigma=Sigma))
}
