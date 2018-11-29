# Description
# Add priors to a plot given the input hyperparameters and number of priors to
# add.
#
# Example calls(s)
#
# bayDem_plotSamples(hp,N,y)
#
# Input(s)
#   Name    Type           Description
#   samps   list           List of samples to plot
#   y       vector         Points at which to evaluate curves
#   rPlot   boolean        [optional] TRUE if the plot is the rate and FALSE if
#                          the plot is the probability density
#
# Output(s)
#   Name    Type           Description
#   NA      NA             NA
bayDem_plotSamples <- function(samps, y, rPlot = F, add = F, lineCol = "black", lineWid = 1) {
  N <- length(samps)
  for (n in 1:N) {
    th <- bayDem_gaussMixParamList2Vect(samps[[n]])
    f <- bayDem_calcGaussMixPdf(th, y)
    if (rPlot) {
      r <- bayDem_f2r(f, y[2] - y[1])
      if (add) {
        lines(y[1:(length(f) - 1)], r, col = lineCol, lwd = lineWid)
      } else {
        plot(y[1:(length(f) - 1)], r, type = "l", col = lineCol, lwd = lineWid)
      }
    } else {
      if (add) {
        lines(y, f, col = lineCol, lwd = lineWid)
      } else {
        plot(y, f, type = "l", col = lineCol, lwd = lineWid)
      }
    }
  }
}
