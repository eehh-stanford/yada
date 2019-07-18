#' @title Convert theta_y from vector to list format
#'
#' @description \code{theta_yVect2List} converts the variable theta_y formatted as a vector to a list. The ordering of the vector is [rho, a, r, b, tau, s].
#'
#' @param theta_y theta_y as a vector
#' @param hp Hyperparameters (needed to guide conversion)
#'
# @keywords
#' @export
#'
#' @return theta_y as a list (theta_y_list)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

theta_yVect2List <- function(theta_y,hp) {
  
  theta_y_list <- list(paramModel=hp$paramModel)
  corr   <- !grepl('uncorr',theta_y_list$paramModel) # Use correlations?
  hetero <-  grepl('hetero',theta_y_list$paramModel) # Use heteroskedasticity?

  last <- 0 # The last index added

  if(hp$J > 0) {
    number <- hp$J # Number of elements to add

    # Add ordinal variables
    theta_y_list$rho <- theta_y[(last+1):(last+number)]
  }

  if(hp$K > 0) {
    # Add continuous variables
    last <- last + number
    number <- hp$K
    theta_y_list$a <- theta_y[(last+1):(last+number)]

    last <- last + number
    theta_y_list$r <- theta_y[(last+1):(last+number)]

    last <- last + number
    theta_y_list$b <- theta_y[(last+1):(last+number)]
  }

  if(hp$J > 0) {
    # Add tau
    tau <- list()
    for(j in 1:hp$J) {
      last <- last + number
      number <- hp$M[j]
      tau[[j]] <- theta_y[(last+1):(last+number)]
    }
    theta_y_list$tau <- tau
  }

  # Add covariance matrix
  last <- last + number
  if(corr) {
    number <- choose(hp$K+hp$J,2) + hp$K + hp$J
    theta_y_list$Sigma <- s2Sigma(theta_y[(last+1):(last+number)])
  } else {
    number <- hp$K + hp$J
    theta_y_list$Sigma <- diag(theta_y[(last+1):(last+number)])
  }

  # Add covariance matrix scale parameter (gamma)
  if(hetero) {
    last <- last + number
    theta_y_list$gamma <- theta_y[last+1]
  }

  return(theta_y_list)
}
