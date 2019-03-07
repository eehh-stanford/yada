#' @title Adaptive simulated annealing optimization
#'
#' @description This function implements adaptive simulated annealing, in which the adaptation is to modify the proposal distribution using the approach in Harrio, Saksman, and Tamminen (2001), An adaptive Metropolis algorithm, in Bernoulli, Volume 7, Number 2, pages 223 through 242. 
#'
#' @param costFunc The cost function to be minimized
#' @param X_0 The starting point for sampling
#' @param ... Further arguments to be passed to costFunc
#' @param control A list of control parameters. See Details.
#'
#' @return A list containing the samples along with summary information
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
adaptiveSA <- function(costFunc,X_0,...,control=list()) {

  # If necessary, set control defaults
  if(!('sampControls' %in% names(control))) {
    sampControls <- list()
  } else {
    sampControls <- control$sampControls
  }

  if(!('T0' %in% names(control))) {
    control$T0 <- 10
  }

  if(!('tempRed' %in% names(control))) {
    control$tempRed <- .95
  }

  if(!('costTol' %in% names(control))) {
    control$costTol <- 1e-7
  }

  if(!('saveFile' %in% names(control))) {
    doSave <- F
  } else {
    doSave <- T
  }

  sampList <- list()
  # First temperature
  sampControls <- control$sampControls
  temp <- control$T0
  sampList[[1]] <- saMetrop(costFunc,X_0,temp,...,control=sampControls)
  if(doSave) {
    saveRDS(sampList,control$saveFile)
  }

  done <- F
  n <- 2
  while(!done) {
    temp <- temp * control$tempRed
    numSamp <- nrow(sampList[[n-1]]$X_mat)
    X_0 <- sampList[[n-1]]$X_mat[,numSamp]
    sampControls$C <- cov(t(sampList[[n-1]]$X_mat)) * control$tempRed
    sampList[[n]] <- saMetrop(costFunc,X_0,temp,...,control=sampControls)
    if(doSave) {
      saveRDS(sampList,control$saveFile)
    }
    done <- sd(sampList[[n]]$cost) <= control$costTol
    n <- n + 1
  }
  return(list(sampList=sampList,control=control))
}
