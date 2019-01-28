#' @title Save a sample of theta_y from the adaptive Metropolis algorithm
#'
#' @description \code{save_theta_y_sample} saves to file the result of sampling theta_y using the Metropolis algorithm
#'
# @details
#'
#' @param fileName The save file
#' @param thetaList A list of samples of theta_y (samples are a list, not a vector)
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param hp Hyperparameters
#' @param varNames Vector of variable names [optional]
#' @param known List object of known simulation values [optional]
#'
# @keywords
#' @export
#'
# @examples
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
# @references

#' @export
save_theta_y_sample <- function(fileName,thetaList,logLikVect,x,Y,hp,varNames=NA,known=NA) {
  bundle <- bundle_theta_y_sample(thetaList,logLikVect,x,Y,hp,varNames,known)
  saveRDS(bundle,fileName)
}
