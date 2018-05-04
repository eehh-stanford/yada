# Description
#   Calculate the log-likelihood for a given set of radiocarbon measurements
#   a parameterization, th, of the (truncated) Gaussian mixture. This can be
#   generalized in the future beyond Gaussian mixtures. The measurement matrix
#   (see bayDem_calcMeasMatrix) is input do reduce computation.
#
# Example calls(s)
#
#   logLik <_ bayDem_calcLogLik(ygrid,th,M)
# 
# Input(s)
#   Name          Type           Description
#   ygrid         vector         The locations of the numeric integration over
#                                calendar dates
#   th            vector-like    The parameterization of the population pdf
#                                (see bayDem_sampleGaussMix)
#   M             matrix         [numMeas x numGrid] The measurement matrix (see
#                                bayDem_calcMeasMatrix)
#
# Output(s)
#   Name          Type           Description
#   logLik        scalar         The log-likelihood, log( p(D|th,alpha) )

bayDem_calcLogLik <- function(ygrid,th,M) {
    fth <- bayDem_calcGaussMixPdf(th,ygrid)
    likVect <- M %*% fth
    logLik <- sum(log(likVect))
    return(logLik)
}
