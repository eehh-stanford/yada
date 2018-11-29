# Description
#   For the input pdf matrix fMat, with dimensions numSamp x numGrid, calculate
#   the "quantiles" at each grid point. The default quantiles are
#   c(.025,.5,.975), yielding the 2.5%, 50% , and 97.5% values.
#
# Example calls(s)
#
#   Q <- bayDem_calcQuantiles(fMat)
#
# Input(s)
#   Name    Type           Description
#   fMat    matrix         The pdf matrix [numSamp x numGrid]
#   probs   vector         [optional] Vector of probabilities for quantile
#                          calculation
#
# Output(s)
#   Name    Type           Description
#   Q       matrix         Matrix of quantiles [length(probs) x numGrid]
bayDem_calcQuantiles <- function(fMat, probs = c(.025, .5, .975)) {
  numQuant <- length(probs)
  numSamp <- dim(fMat)[1]
  numGrid <- dim(fMat)[2]
  Q <- matrix(NA, numQuant, numGrid)
  for (g in 1:numGrid) {
    Q[, g] <- quantile(fMat[, g], probs = probs)
  }
  return(Q)
}
