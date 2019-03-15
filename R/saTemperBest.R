#' @title Get best value from call to saTemper
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
saTemperBest <- function(temper) {
  bestChain <- which.min(unlist(lapply(temper$chainList,function(x){min(x$costVect)})))
  theta <- temper$chainList[[bestChain]]$X_mat[,which.min(temper$chainList[[bestChain]]$costVect)]
  return(theta)
}
