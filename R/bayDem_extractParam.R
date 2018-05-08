# Description
#   Calculate the probability density function (pdf) for a (truncated)
#   Gaussian mixture with K components:
#
#   f = [p(y_1|th)]
#       [p(y_2|th)]
#       [   ...   ]
#       [   ...   ]
#
# Example calls(s)
#
#   TH <- bayDem_extractParam(fit,hp)
#   samps <- bayDem_extractParam(fit,hp,TRUE)
# 
# Input(s)
#   Name    Type           Description
#   fit     rstan          The rstan fit object
#   hp      list           The hyperparameters
#   asList  logical        [default: FALSE] Whether to return a list of lists
#                          (stan formatting) or a matrix (more efficient for
#                          calculations)
#
# Output(s)
#   Name    Type           Description
#   TH      matrix         [numSamp x numParam] Matrix of parameters from the
#                          posterior sample in fit (if asList is FALSE)
#   samps   list           List of parameters from the posterior sample in fit
#                          (if asList is TRUE)
bayDem_extractParam <- function(fit,hp,asList=F) {
  X <- as.data.frame(fit)
  numSamp <- dim(X)[1]
   
  if(hp$fitType=='gaussmix') {
    numParam <- 2 + 3*hp$K
    piVect <- rep(NA,hp$K)
    muVect <- rep(NA,hp$K)
    sigVect <- rep(NA,hp$K)
    TH <- rep(NA,numSamp,numParam)
    varList <- c()
    for(k in 1:hp$K) {
      varList[k] <- paste('pi[',as.character(k),']',sep='')
      varList[k + hp$K] <- paste('mu[',as.character(k),']',sep='')
      varList[k + 2*hp$K] <- paste('sig[',as.character(k),']',sep='')
    }
    TH <- as.matrix(X[,varList])
    TH <- cbind(rep(hp$ymin,numSamp),rep(hp$ymax,numSamp),TH)
    colnames(TH) <- c('ymin','ymax',varList)

    if(!asList) {
      class(TH) <- hp$fitType # Set the class to aid subsequent evaluation
      return(TH)
    } else {
      samps <- list()
      for(ss in 1:numSamp) {
        samps[[ss]] <- bayDem_gaussMixParamVect2List(TH[ss,])
      }
      return(samps)
    }
  } else {
    stop(paste('Unrecognized fit type:',hp$fitType))
  }
}
