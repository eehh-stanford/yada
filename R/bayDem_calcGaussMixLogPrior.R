# Description
#   For the (truncated) Gaussian mixture, calculate the (log) prior probability
#   of the parameter theta given the hyperparameters. Since the mean is
#   un-informative (it is drawn from the uniform distribution on ymin to ymax)
#   this component of the prior probability is not calculated by default. This
#   has no impact on the Metropolis-Hastings sampling.
#
# Example calls(s)
#   logP <- bayDem_calcGaussMixLogPrior(th,gmHyper)
# 
# Input(s)
#   Name       Type     Description
#   th         vector   The parameter vector (see bayDem_calcGaussMixPdf)
#   gmHyper    list     The hyperparameters:
#                         $z$zAlpha      -- Parameter for the Dirichlet draw
#                         $sig$sigAlpha  -- Shape parameter for the gamma draw
#                         $sig$sigBeta   -- Rate parameter for the gamma draw
#                         $sig$sigOffset -- Offset for the gamma draw
#
# Output(s)
#   Name       Type     Description
#   logP       scalar   The (log) prior probability, log p(theta|alpha)

bayDem_calcGaussMixLogPrior <- function(th,gmHyper) {
    J <- (length(th)-2)/3 # Number of  mixtures

    # Sample for mixture proportions (Z is [numSamp x 2]
    pLog <- log(MCMCpack::ddirichlet(th[3:(2+J)],rep(gmHyper$z$zAlpha,J)))

    pLog <- pLog + sum(dgamma(th[(3+J):(2+2*J)]-gmHyper$sig$sigOffset,shape=gmHyper$sig$sigAlpha,rate=gmHyper$sig$sigBeta))
    return(pLog)
}
