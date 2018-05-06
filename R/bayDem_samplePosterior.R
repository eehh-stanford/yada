# Description
# Get a random sample of the posteriors from the stan fit.
#
# Example calls(s)
#
#   samp <- bayDem_samplePosterior(fit,N)
# 
# Input(s)
#   Name    Type           Description
#   fit     rstan          The rstan fit object
#   hp      list           The hyperparameters
#   N       integer        Number of samples
#
# Output(s)
#   Name    Type           Description
#   samps   list           N samples from the posterior for the given fit type.
#                          samps is a list of lists, which is the format
#                          expected by stan for initializing chains. Details for
#                          each fit type:
#
#                gaussmix  sig [K x 1]  -- vector of standard deviations
#                          mu  [K x 1]  -- vector of means
#                          pi  [K x 1]  -- vector of mixture proportions
#                          ymin         -- minimum calendar date
#                          ymax         -- maximum calendar date

bayDem_samplePosterior <- function(fit,hp,N) {
  numChains <- fit@sim$chains
  sampsPerChain <- fit@sim$iter
  numTot <- numChains * sampsPerChain
  ind <- sample.int(numTot,size=N)
  samps <- list()
  for(n in 1:N) {
    cc <- 1 + floor((ind[n]-1) / sampsPerChain) # The chain
    ss <- ((ind[n]-1) %% sampsPerChain) + 1 # Sample within chain
    if(hp$fitType == 'gaussmix') {
      piVect <- rep(NA,K)
      muVect <- rep(NA,K)
      sigVect <- rep(NA,K)
      for(k in 1:hp$K) {
        piVect[k]  <- fit@sim$samples[[cc]][[paste('pi[',as.character(k),']',sep='')]][ss]
        muVect[k]  <- fit@sim$samples[[cc]][[paste('mu[',as.character(k),']',sep='')]][ss]
        sigVect[k] <- fit@sim$samples[[cc]][[paste('sig[',as.character(k),']',sep='')]][ss]
      }
      samp <- list(fitType=hp$fitType,ymin=hp$ymin,ymax=hp$ymax,pi=piVect,mu=muVect,sig=sigVect)
    } else {
      stop(paste('Unrecognized fit type:',hp$fitType))
    }
    samps[[n]] <- samp
  }
  return(samps)
}

