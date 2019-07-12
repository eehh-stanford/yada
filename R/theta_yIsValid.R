#' @title Check if theta_y is valid for sampling
#'
#' @description \code{theta_yIsValid} checks if the input theta_y (as a list) is valid for sampling. In particular, ensure that (1) tau is ordered, (2) all correlations are between -1 and 1, (3) Sigma is positive definite, (4) and [optional] the likelihood is finite. 
#'
#' @param theta_y_list theta_y as a list
#' @param hp Hyperparameters
#' @param forceFinite Whether to ensure likelihood is finite [optional]
#' @param x Vector of indepenent variable observations [optional]
#' @param Y Matrix of dependent variable observations [optional]
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return Boolean indicating whether theta_y is valid
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

theta_yIsValid <- function(theta_y_list,hp,forceFinite=F,x=x,Y=Y) {
  # Check if gamma (covariance scaling) is between 0 and 1
  if(theta_y_list$gamma <= 0 || theta_y_list$gamma > 1) {
    return(F)
  }

  if(hp$J > 0) {
    # Check if tau is ordered
    for(j in 1:hp$J) {
      if(is.unsorted(theta_y_list$tau[[j]])) {
        return(F)
      }
    }
  }

  # Check if Sigma is positive definite
  if(!matrixcalc::is.positive.definite(theta_y_list$Sigma)) {
    return(F)
  }

  # If necessary, ensure likelihood is finite
  if(forceFinite) {
   ll <- calcLogLik_theta_y(theta_y_list,x,Y,hp)
   if(is.finite(ll)) {
     return(T)
   } else {
     return(F)
   }
  } else {
    return(T)
  }
}
