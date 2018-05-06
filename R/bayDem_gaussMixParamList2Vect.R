# Description
#   For the input, parameterized (truncated) Guassian convert from a list
#   representation (required by stan) to a vector representation (which is
#   more efficient for calculations and the form expected by, e.g.,
#   bayDem_calcGaussMixPdf.R).
#
# Example calls(s)
#
#   th <- bayDem_gaussMixParamList2Vect(samp)
# 
# Input(s)
#   Name    Type           Description
#   samp    list           A list that parameterizes the (truncated) Gaussian
#                          mixture. See bayDem_samplePrior.R.
#
# Output(s)
#   Name    Type           Description
#   th      vector         A vector that parameterizes the (truncated) Gaussian
#                          mixture. See bayDem_calcGaussMixPdf.R
bayDem_gaussMixParamList2Vect <- function(samp) {
  K <- length(samp$sig)
  th <- c(samp$ymin,samp$ymax,samp$pi,samp$mu,samp$sig)
  return(th)
}
