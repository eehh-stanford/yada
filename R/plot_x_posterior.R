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
  plot(1, type="n", xlab="x", ylab="density", xlim=c(min(post$x),max(post$x)),ylim=c(0,max(post$density)))
  xv <- c(post$x[(post$x < post$xlo)],post$xlo)
  fv <- c(post$density[(post$x < post$xlo)],post$flo)
  polygon(c(xv,rev(xv)),c(rep(0,length(fv)),rev(fv)),col='grey',border=NA)
  xv <- c(post$xhi,post$x[(post$x > post$xhi)])
  fv <- c(post$fhi,post$density[(post$x > post$xhi)])
  polygon(c(xv,rev(xv)),c(rep(0,length(fv)),rev(fv)),col='grey',border=NA)
  lines(c(1,1)*post$xmean,c(0,post$fmean),lty=2,lwd=3,col='grey')
  lines(c(1,1)*post$xmed,c(0,post$fmed),lwd=3,col='grey')
  if(!is.null(post$xknown)) {
    lines(c(1,1)*post$xknown,c(0,post$fknown),lwd=3,col='black')
  }
  lines(post$x,post$density,type='l',lwd=3)
  if(!is.null(post$xknown)) {
    legend('topright',legend=c('Actual','Median','Mean','+/- 2.5%'),lty=c(1,1,2,NA),col=c('black','grey','grey','grey'),pch=c(NA,NA,NA,22),pt.bg=c(NA,NA,NA,'grey'),lwd=c(2,2,2,NA),pt.cex=c(1,1,1,2))
  } else {
    legend('topright',legend=c('Median','Mean','+/- 2.5%'),lty=c(1,2,NA),col='grey',pch=c(NA,NA,22),pt.bg=c(NA,NA,'grey'),lwd=c(2,2,NA),pt.cex=c(1,1,2))
  }
}
