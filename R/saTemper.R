#' @title Parallel tempering for function minimization
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
saTemper <- function(costFunc,init,tempVect=NA,...,control=list()) {

  haveChains <- is.list(init)
  if(!haveChains) {
    numTemp <- length(tempVect)
    X_0 <- init
    X_0 <- t(matrix(rep(X_0,numTemp),ncol=numTemp))
  } else {
    tempVect <- init$tempVect # tempVect is ignored if input
    numTemp <- length(tempVect)
    control0 <- control
    control <- init$control
  }
  # If necessary, set control defaults
  if(!haveChains) {
    if(!('metropControl' %in% names(control))) {
      metropControl <- list()
    } else {
      metropControl <- control$metropControl
    }
  }

  if(haveChains) {
    if(('metropControl' %in% names(control0))) {
      metropControl <- control0$metropControl
    } else {
      metropControl <- init$metropControl
    }
  }

  if(!haveChains && !('numCycles' %in% names(control))) {
    control$numCycles <- 100 # Number of Metropolis-Hastings / swap cycles
  }

  if(haveChains && ('numCycles' %in% names(control0))) {
    control$numCycles <- control0$numCycles
  }

  if(!haveChains && !('verbose' %in% names(control))) {
    control$verbose <- F
  }

  if(haveChains && ('verbose' %in% names(control0))) {
    control$verbose <- control0$verbose
  }


  if(!haveChains) {
    chainList <- list()
  } else {
    chainList <- init$chainList
  }

  # Run each chain once [if not picking up from a previous set of runs]
  if(!haveChains) {
    for(m in 1:numTemp) {
      chainList[[m]] <- saMetrop(costFunc,as.vector(X_0[m,]),tempVect[m],...,control=metropControl)
    }
    # If this is the first run and the number of cycles is 1, return without swapping
    if(control$numCycles == 1) {
      return(list(chainList=chainList,control=control,tempVect=tempVect))
    }
  }


  done <- F
  #swapList <- list()
  if(!haveChains) {
    cycle <- 1
  } else {
    cycle <- 0
  }
  while(!done) {
    # Attempt swaps
    if(control$verbose) {
      print('**--**--**')
      if(!haveChains) {
        print(cycle)
      } else {
        print(cycle + 1)
      }
      print(min(unlist(lapply(chainList,function(x){min(x$costVect)}))))
      state0 <- unlist(lapply(chainList,function(x){x$costVect[length(x$costVect)]}))
    }
    for(m in 1:(numTemp-1)) {
      c_i <- chainList[[m]]$costVect[length(chainList[[m]]$costVect)]
      c_j <- chainList[[m+1]]$costVect[length(chainList[[m+1]]$costVect)]
      T_i <- tempVect[m]
      T_j <- tempVect[m+1]
      alpha_swap <- exp(-(c_j-c_i)*(1/T_i-1/T_j))
      if(alpha_swap >= 1) {
        accept <- T
      } else {
        accept <- runif(1) <= alpha_swap
      }
      #if(length(swapList) == m-1) {
      #  swapList[[m]] <- accept
      #} else {
      #  swapList[[m]] <- c(swapList[[m]],accept)
      #}

      if(accept) {
        X_i <- chainList[[m]]$X_mat[,ncol(chainList[[m]]$X_mat)]
        X_j <- chainList[[m+1]]$X_mat[,ncol(chainList[[m+1]]$X_mat)]
        c_i <- chainList[[m]]$costVect[length(chainList[[m]]$costVect)]
        c_j <- chainList[[m+1]]$costVect[length(chainList[[m+1]]$costVect)]
        a_i <- chainList[[m]]$acceptVect[length(chainList[[m]]$acceptVect)]
        a_j <- chainList[[m+1]]$acceptVect[length(chainList[[m+1]]$acceptVect)]

        chainList[[m]]$X_mat[,ncol(chainList[[m]]$X_mat)] <- X_j
        chainList[[m+1]]$X_mat[,ncol(chainList[[m+1]]$X_mat)] <- X_i
        chainList[[m]]$costVect[length(chainList[[m]]$costVect)] <- c_j
        chainList[[m+1]]$costVect[length(chainList[[m+1]]$costVect)] <- c_i
        chainList[[m]]$acceptVect[length(chainList[[m]]$acceptVect)] <- a_j
        chainList[[m+1]]$acceptVect[length(chainList[[m+1]]$acceptVect)] <- a_i
      }
      #if(control$verbose) {
      #  print('--')
      #  print(m)
      #  print(alpha_swap)
      #  print(accept)
      #}
    }
    if(control$verbose) {
      state1 <- unlist(lapply(chainList,function(x){x$costVect[length(x$costVect)]}))
      print(rbind(state0,state1))
    }

    for(m in 1:numTemp) {
      chainList[[m]] <- saMetrop(costFunc,chainList[[m]],tempVect[m],...,control=metropControl)
    }
    cycle <- cycle + 1
    done <- cycle == control$numCycles
  }
  return(list(chainList=chainList,control=control,tempVect=tempVect))
  #return(list(chainList=chainList,swapList=swapList,control=control,tempVect=tempVect))
}
