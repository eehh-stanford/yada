#' @title Bundle result of Metropolis sampling into a single convenient list
#'
#' @description \code{bundles_theta_y_sample} bundles theta_yList, x, Y, hp, varNames, and known into a single convenient list. This bundled result is used for displaying and saving results.
#'
# @details
#'
#' @param theta_yList A list of samples of theta_y (samples are a list, not a vector)
#' @param logLikeVect A vector of likelihoods calculated for each sample
#' @param x Vector of indepenent variable observations
#' @param Y Matrix of dependent variable observations
#' @param hp Hyperparameters
#' @param varNames Vector of variable names [optional]
#' @param known List object of known simulation values [optional]
#'
# @keywords
#' @export
#'
#' @return The bundled list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
# @references

#' @export
bundle_theta_y_sample <- function(theta_yList,logLikVect,x,Y,hp,varNames=NA,known=NA) {
  bundle <- list(theta_yList=theta_yList,logLikVect=logLikVect,Y=Y,x=x,hp=hp)
  if(!all(is.na(known))) {
    bundle$known <- known
    bundle$haveKnown <- T
  } else {
    bundle$haveKnown <- F
  }

  if(!all(is.na(varNames))) {
    bundle$varNames <- varNames
    bundle$haveNames <- T
  } else {
    bundle$haveNames <- F
  }

  return(bundle)
}
