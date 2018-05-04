# Description
#   Calculate the log-likelihood for a given set of radiocarbon measurements
#   given the parametrized wavelet specified by the augData, the augmented
#   data.
#
# Example calls(s)
#
# 
# Input(s)
#   Name          Type           Description
#   M             Matrix         [Nmeas x Ngrid] The measurement matrix (see
#                                bayDem_calcMeasMatrix)
#   augData       List           The augmented data specifying the 
#
# Output(s)
#   Name          Type           Description
#   logLik        scalar         The log-likelihood, log( p(D|th,alpha) )

bayDem_calcWaveletLogLik <- function(M,augData,dy=1) {
    f <- bayDem_calcWaveletPdf(augData,dy)
    f <- as.matrix(f,n.col=1)
    likVect <- M %*% f
    logLik <- sum(log(likVect))
    return(logLik)
}
