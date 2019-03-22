#' @title Plot the result of a call to analyze_x_posterior
#'
#' @description \code{plot_x_posterior} plots the result of a call to analyze_x_posterior. In particular, the density, median, mean, and 2.5% confidence interval locations are shown, as is the actual age (if known).
#'
# @details
#'
#' @param post The output of analyze_x_posterior
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
plot_x_posterior <- function(post) {
  locut <- post$xlolo
  if(locut <=1) {
    locut <- 0
  }
  hicut <- post$xhihi
  x <- post$x
  ind <- locut <= x & x <= hicut
  x <- x[ind]
  f <- post$density[ind]

  plot(1, type="n", xlab="x", ylab="density", xlim=c(min(locut),max(hicut)),ylim=c(0,max(f)))
  xv <- c(x[(x < post$xlo)],post$xlo)
  fv <- c(f[(x < post$xlo)],post$flo)
  polygon(c(xv,rev(xv)),c(rep(0,length(fv)),rev(fv)),col='grey',border=NA)
  xv <- c(post$xhi,x[(x > post$xhi)])
  fv <- c(post$fhi,f[(x > post$xhi)])
  polygon(c(xv,rev(xv)),c(rep(0,length(fv)),rev(fv)),col='grey',border=NA)
  lines(c(1,1)*post$xmean,c(0,post$fmean),lty=2,lwd=3,col='grey')
  if(!is.null(post$xknown)) {
    lines(c(1,1)*post$xknown,c(0,post$fknown),lwd=3,col='black')
  }
  lines(x,f,type='l',lwd=3)
  if(!is.null(post$xknown)) {
    legend('topright',legend=c('Actual','Mean','+/- 2.5%'),lty=c(1,1,2,NA),col=c('black','grey','grey'),pch=c(NA,NA,22),pt.bg=c(NA,NA,'grey'),lwd=c(2,2,NA),pt.cex=c(1,1,2))
  } else {
    legend('topright',legend=c('Mean','+/- 2.5%'),lty=c(2,NA),col='grey',pch=c(NA,22),pt.bg=c(NA,'grey'),lwd=c(2,NA),pt.cex=c(1,2))
  }
}
