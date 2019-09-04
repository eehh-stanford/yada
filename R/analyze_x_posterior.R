#' @title Analyze an input density by calculating the mean and confidence intervals
#'
#' @description \code{analyze_x_posterior} analyzes an input density given by the input pair of vectors xv and fv by calculating the mean and confidence intervals for a handful of quantiles (0.001, 0.025, 0.5, 0.975, and 0.999). xv is assumed to be evenly spaced. Optionally, a known age (xknown) can be input, for which the density is estimated.
#'
# @details
#'
#' @param xv A vector of evenly spaced ages at the density, fv, is calculated
#' @param fv A vector of densities
#' @param xknown Known age [optional]
#'
#'
#' @return A list containing the results of the analysis
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
analyze_x_posterior <- function(xv,fv,xknown=NA) {
  if(xv[1] != 0) {
    stop('First element of xv should be zero')
  }
  dx <- xv[2] - xv[1]
  # fv is the same as p_x
  fv <- fv / sum(fv) / dx
  Fv <- c(0,dx*cumsum(fv)[1:(length(xv)-1)])
  
  q <- 0.001
  n <- max(which(Fv <= q))
  xlolo <- xv[n] + dx*(q-Fv[n])/(Fv[n+1]-Fv[n])
  flolo <- fv[n] + (xlolo-xv[n])*(fv[n+1]-fv[n])/dx

  q <- 0.025
  n <- max(which(Fv <= q))
  xlo <- xv[n] + dx*(q-Fv[n])/(Fv[n+1]-Fv[n])
  flo <- fv[n] + (xlo-xv[n])*(fv[n+1]-fv[n])/dx

  q <- 0.5
  n <- max(which(Fv <= q))
  xmed <- xv[n] + dx*(q-Fv[n])/(Fv[n+1]-Fv[n])
  fmed <- fv[n] + (xmed-xv[n])*(fv[n+1]-fv[n])/dx

  q <- 0.975
  n <- max(which(Fv <= q))
  xhi <- xv[n] + dx*(q-Fv[n])/(Fv[n+1]-Fv[n])
  fhi <- fv[n] + (xhi-xv[n])*(fv[n+1]-fv[n])/dx

  q <- 0.999
  n <- max(which(Fv <= q))
  xhihi <- xv[n] + dx*(q-Fv[n])/(Fv[n+1]-Fv[n])
  fhihi <- fv[n] + (xhihi-xv[n])*(fv[n+1]-fv[n])/dx

  xmean <- sum(fv*xv)*dx
  n <- max(which(xv < xmean))
  fmean <- fv[n] + (xmean-xv[n])*(fv[n+1]-fv[n])/dx

  returnList <- list(x=xv,density=fv,dx=dx,xlolo=xlolo,xlo=xlo,xmed=xmed,xhi=xhi,xhihi=xhihi,xmean=xmean,flo=flo,fmed=fmed,fhi=fhi,fmean=fmean)
  if(!is.na(xknown)) {
    # Calculate the expectation (over the density) of (x-xknown)^2
    expectedSqrErr <- sum(dx*fv * (xv-xknown)^2)
    if(xknown == 0) {
      fknown <- fv[1]
    } else {
      n <- max(which(xv < xknown))
      fknown <- fv[n] + (fv[n+1] - fv[n]) * (xknown-xv[n]) / dx
      Fknown <- Fv[n] + (Fv[n+1] - Fv[n]) * (xknown-xv[n]) / dx
    }
    returnList$expectedSqrErr <- expectedSqrErr
    returnList$xknown <- xknown
    returnList$fknown <- fknown
    returnList$Fknown <- Fknown
  }
  return(returnList)
}
