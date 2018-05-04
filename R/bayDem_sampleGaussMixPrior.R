# Description
#   Sample from the prior distribution p(th|alpha) for the (truncated) Gaussian
#   mixture (see also bayDem_calcGaussMixPriorParam).
#
# Example calls(s)
#
#   TH <- bayDem_sampleGaussMixPrior(gmPriorParam,numSamp)
# 
# Input(s)
#   Name            Type     Description
#   gmPriorParam    list     A list that parametrizes the model parameter
#                            theta, with the following components:
#                            ymin -- Mininum calendar date
#                            ymax -- Maximum calendar date
#                            startOffset -- constant added to gamma draw
#                            modeOffset  -- location of mode for gamma draw
#                            gammaAlpha  -- alpha parameter for gamma draw
#                            gammaBeta   -- beta parameter for gamma draw
#                            dirichParam -- dirichlet parameters for z
#   numSamp         int      Number of samples
#
# Output(s)
#   Name            Type     Description
#   TH              matrix   A matrix of samples [6 x numSamp] with the
#                            following rows:
#                            z1   -- Mixture proportion of first mixture
#                            z2   -- Mixture proportion of second mixture
#                                    (z2 = 1 - z1)
#                            mu1  -- Mean of first normal mixture
#                            sig1 -- Standard deviation of first normal mixture
#                            mu2  -- Mean of second normal mixture
#                            sig2 -- Standard deviation of second normal mixture
#                            ymin -- Left truncation point
#                            ymax -- Right truncation point

bayDem_sampleGaussMixPrior <- function(gmPriorParam,numSamp) {
    # Sample for mixture proportions (Z is [numSamp x 2]
    Z <- MCMCpack::rdirichlet(numSamp,c(gmPriorParam$dirichParam,gmPriorParam$dirichParam))
    # Sample for means (MU is [2*numSamp x 1])
    MU <- runif(2*numSamp,gmPriorParam$ymin,gmPriorParam$ymax)

    # Sample for standard deviations (SIG is [2*numSamp x 1])
    # Samples for standards deviations
    SIG <- rgamma(2*numSamp,shape=gmPriorParam$gammaAlpha,rate=gmPriorParam$gammaBeta) + gmPriorParam$startOffset

    # Truncation limits
    yMinVect <- rep(gmPriorParam$ymin,numSamp)
    yMaxVect <- rep(gmPriorParam$ymax,numSamp)

    # Create the [8 x numSamp] output matrix, TH
    TH <- rbind(Z[,1],Z[,2],MU[1:numSamp],SIG[1:numSamp],MU[(numSamp+1):(2*numSamp)],SIG[(numSamp+1):(2*numSamp)],yMinVect,yMaxVect)
    rownames(TH) <- c('z1','z2','mu1','sig1','mu2','sig2','ymin','ymax')
    return(TH)
}
