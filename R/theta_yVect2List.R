#' @title Convert theta_y from vector to list format
#'
#' @description \code{theta_yVect2List} converts the variable theta_y formatted as a vector to a list. The ordering of the vector is [alpha, rho, a, r, b, tau, s].
#'
#' @param theta_y theta_y as a vector
#' @param hp Hyperparameters (needed to guide conversion)
#' @param theta_y_list theta_y as a list
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return theta_y as a vector
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

theta_yVect2List <- function(theta_y,hp) {
  theta_y_list <- list(paramModel=hp$paramModel)
  
  last <- 0 # The last index added
  number <- hp$J # Number of elements to add

  # Add ordinal variables
  theta_y_list$alpha <- theta_y[(last+1):(last+number)]

  last <- last + number
  theta_y_list$rho <- theta_y[(last+1):(last+number)]

  # Add continuous variables
  last <- last + number
  number <- hp$K
  theta_y_list$a <- theta_y[(last+1):(last+number)]

  last <- last + number
  theta_y_list$r <- theta_y[(last+1):(last+number)]

  last <- last + number
  theta_y_list$b <- theta_y[(last+1):(last+number)]

  # Add tau
  tau <- list()
  for(j in 1:hp$J) {
    last <- last + number
    number <- hp$M[j]
    tau[[j]] <- theta_y[(last+1):(last+number)]
  }
  theta_y_list$tau <- tau

  # Add covariance matrix
  last <- last + number
  number <- choose(hp$K+hp$J,2)
  theta_y_list$Sigma <- corrVect2covMat(theta_y[(last+1):(last+number)])

  return(theta_y_list)
}
