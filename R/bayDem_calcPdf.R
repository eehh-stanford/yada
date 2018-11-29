# @keywords
#' @export
# Description
#   Calculate the probability density function (pdf) for the input
#   vector-like parameterization. The class of the pdf is either input or
#   determined from the input.
#
# Example calls(s)
#
#   f <- bayDem_calcPdf(th,y)
#   f <- bayDem_calcPdf(th,y,'gaussmix')
#
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vector-like object that parameterizes the pdf
#   y       vector         A vector of grid points at which to evaluate the pdf
#   fitType string         The parameterization for this pdf. If it is not
#                          input, it assumed to be the class of th.
#
# Output(s)
#   Name    Type           Description
#   f       vector         Output pdf
bayDem_calcPdf <- function(th, y, fitType = NA) {
  if (is.na(fitType)) {
    fitType <- class(th)
  }

  if (fitType == "gaussmix") {
    return(bayDem_calcGaussMixPdf(th, y))
  } else {
    stop(paste("Unrecognized fit type:", hp$fitType))
  }
}
