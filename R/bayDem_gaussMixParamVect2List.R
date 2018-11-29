# Description
#   For the input, parameterized Guassian mixture convert from a vector
#   representation (which is more efficient for calculations and the form
#   expected by, e.g., bayDem_calcGaussMixPdf.R) to a vector representation
#   (required by stan).
#
# Example calls(s)
#
#   samp <- bayDem_gaussMixParamVect2List(th)
#
# Input(s)
#   Name    Type           Description
#   th      vector         A vector that parameterizes the Gaussian mixture.
#                          See bayDem_calcGaussMixPdf.R
#
# Output(s)
#   Name    Type           Description
#   samp    list           A list that parameterizes the Gaussian mixture. See
#                          bayDem_samplePrior.R.
bayDem_gaussMixParamVect2List <- function(th) {
  K <- length(th) / 3 # Number of  mixtures
  samp <- list(pi = th[1:K], mu = th[(1 + K):(2 * K)], sig = th[(1 + 2 * K):(3 * K)])
  return(samp)
}
