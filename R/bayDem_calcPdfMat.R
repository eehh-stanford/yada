# Description
#   Calculate the probability density function (pdf) for the input matrix,
#   which has dimensions numSamp x numParam, where numSamp is the number of
#   samples and numParam is the number of parameters for the fit type. The
#   output is a pdf matrix with dimensions numSamp x numGrid, where numGrid is
#   the length of the input vector y. The class of the pdf (fit type) is either
#   input or determined from the input.
#
# Example calls(s)
#
#   fMat <- bayDem_calcPdfMat(TH,y)
#   fMat <- bayDem_calcPdfMat(TH,y,'gaussmix')
# 
# Input(s)
#   Name    Type           Description
#   TH      matrix         A matrix of parameters [numSamp x numParam]
#   y       vector         A vector of grid points at which to evaluate the pdf
#                          [numGrid = length(y)]
#   fitType string         The parameterization for this pdf. If it is not
#                          input, it assumed to be the class of TH.
#
# Output(s)
#   Name    Type           Description
#   fMat    matrix         Output pdf matrix [numSamp x numGrid]
bayDem_calcPdfMat <- function(TH,y,fitType=NA) {
  if(is.na(fitType)) {
    fitType <- class(TH)
  }
  numSamp <- dim(TH)[1]
  numGrid <- length(y)
  fMat <- matrix(NA,numSamp,numGrid)

  # Even though bayDem_calcPdfMat merely calls bayDem_calcPdf, generation of
  # the pdf matrix is centralized here so that, in the future, it can be
  # calculated in C or via parallelization.
  for(n in 1:numSamp) {
    fMat[n,] <- bayDem_calcPdf(TH[n,],y,fitType)
  }
  return(fMat)
}
