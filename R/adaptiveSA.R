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

  if(!('tempVect' %in% names(control))) {
    control$tempVect <- 16*(15/16)^(0:60)
    #control$tempVect <- c(2,1.75,1.5,1.25,1,.75,.25)
    #control$tempVect <- c(4,2,1)
  }

  sampList <- list()
  print('xxxx')
  print(1)
  sampControls$temp <- control$tempVect[1]
  sampList[[1]] <- saMetrop(costFunc,X_0,control$tempVect[1],...,control=sampControls)
  saveRDS(sampList,'sampList.rds')
  sampControls 
  for(n in 2:length(control$tempVect)) {
    print('xxxx')
    print(n)
    numSamp <- nrow(sampList[[n-1]]$X_mat)
    X_0 <- sampList[[n-1]]$X_mat[,numSamp]
    sampControls$C <- cov(t(sampList[[n-1]]$X_mat)) * control$tempVect[n] / control$tempVect[n-1]
    sampControls$temp <- control$tempVect[n]
    sampList[[n]] <- saMetrop(costFunc,X_0,control$tempVect[n],...,control=sampControls)
    saveRDS(sampList,'sampList.rds')
  }
  return(list(sampList=sampList,control=control))
}
