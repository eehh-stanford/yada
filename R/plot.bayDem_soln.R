# Description
#   Plot the result of Bayesian inference in the input soln object, which is of
#   class bayDem_soln.
#
# Example calls(s)
#
#   plot(soln)
# 
# Input(s)
#   Name    Type           Description
#   soln    list           The solution, a list-like object of class
#                          bayDem_soln (see bayDem_doInference)
#   y       vector         (optional) The calendar dates at which to evaluate
#                          densities. If y is not input, y is built from the
#                          hyperparameters.
#
# Output(s)
#   Name    Type           Description
#   NA      NA             NA
plot.bayDem_soln <- function(soln,y=NA,th_sim=NA) {
  if(all(is.na(y))) {
    y <- seq(soln$prob$hp$ymin,soln$prob$hp$ymax,by=soln$prob$hp$dy)
  }
  plotSim <- !all(is.na(th_sim))

  TH <- bayDem_extractParam(soln$fit,soln$prob$hp)
  fMat <- bayDem_calcPdfMat(TH,y)
  Q <- bayDem_calcQuantiles(fMat)
  M <- bayDem_calcMeasMatrix(y,soln$prob$phi_m,soln$prob$sig_m,soln$prob$calibDf)
  f_spdf <- colSums(M)
  f_spdf <- f_spdf / sum(f_spdf) / (y[2]-y[1]) # This assumes y is evenly spaced
  plotMaxF <- max(max(Q),f_spdf)
  if(plotSim) {
    f_sim <- bayDem_calcPdf(th_sim,y,soln$prob$hp$fitType)
    plotMaxF <- max(plotMaxF,f_sim)
  }
  plot(y,f_spdf,type='l',col='black',lwd=3,ylim=c(0,plotMaxF),xlab='Calendar Date [AD]',ylab='Density')
  lines(y,Q[1,],col='red',lwd=3,lty=2)
  lines(y,Q[2,],col='red',lwd=3,lty=1)
  lines(y,Q[3,],col='red',lwd=3,lty=2)
  if(plotSim) {
    lines(y,f_sim,col='blue',lwd=3,lty=1)
    legend('topleft',c('Target','Summed','Posterior','Posterior 95%'),lwd=3,lty=c(1,1,1,3),col=c('blue','black','red','red'))
  } else {
    legend('topleft',c('Summed','Posterior','Posterior 95%'),lwd=3,lty=c(1,1,3),col=c('black','red','red'))
  }
}
