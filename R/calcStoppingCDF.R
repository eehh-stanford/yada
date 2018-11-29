#' @title Calculate the empirical cumulative density function for a stopping distribution
#'
#' @description Calculate the empirical cumulative density function (CDF) of a sample of points (y) for use in the discontinuous Kolmogorov-Smirnov test at the points in y.calc. The samples are assumed to all be on the interval 0 to 1 (as for the Martingale Survival Residual (MSR)), and the CDF is assumed to be continuous up to a stopping time at which the transition is guaranteed to occur if it has not occurred already. Whether the stopping time was reached is indicated by the boolean vector halt. Duplicate values of y are acceptable. The optional input "where" indicates whether the CDF should be calculated after the  transition (default), instantaneously before the transition transition, or both; this only matters if y.calc contains values that are also in y (i.e, a calculation is done where a transition occurs in y). If y.calc is not specified, it is assumed to be all the values of y.
#'
# @details
#'
#' @param y A vector of samples
#' @param halt A boolean vector indicating whether the stopping time was reached
#' @param y.calc Points at which to calculate the CDF
#' @param where Whether to calculate the CDF before or after the transition (or both)
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return The CDF evaluated at the unique values of y (or y.calc) as a list with the following entries
#'
#' \describe{
#'   \item{\code{y}}{The unique values of y or y.calc}
#'   \item{\code{before}}{The value of the CDF immediately before each y}
#'   \item{\code{after}}{The value of the CDF immediately after each y}
#' }
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcECDFPointPairs, addECDFLinesToPlot
#
#' @references
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}


calcStoppingCDF <- function(y, halt, y.calc = NA, where = "after") {
  if (length(y.calc) == 1) {
    if (all(is.na(y.calc))) {
      y.calc <- sort(unique(y))
    }
  }

  N <- length(y.calc)

  if (where == "after") {
    calcBefore <- F
    calcAfter <- T
    cdfAfter <- rep(NA, N)
  } else if (where == "before") {
    calcBefore <- T
    calcAfter <- F
    cdfBefore <- rep(NA, N)
  } else if (where == "both") {
    calcBefore <- T
    calcAfter <- T
    cdfAfter <- rep(NA, N)
    cdfBefore <- rep(NA, N)
  } else {
    stop("Unrecognized option for input variable \'where\'")
  }

  numSamples <- length(y)
  for (i in 1:N) {
    y0 <- y.calc[i]
    indLesser <- y < y0

    # Calculate the normalized number of transitions at y0 accounting for
    # all transitions before  y0. This includes a component due to
    # regular events and a component due to halting at the stopping
    # time.
    transBefore <- sum(indLesser & !halt) + sum((y0 - y[indLesser & halt]) / (1 - y[indLesser & halt]))

    # If necessary, calculate the CDF accounting for only transitions before y0
    if (calcBefore) {
      cdfBefore[i] <- transBefore / numSamples
    }

    # If necessary, calculate the CDF accounting for transitions up to and including y0
    if (calcAfter) {
      indEqual <- y == y0
      # Calculate the number of transitions exactly at y0. Only regular events
      # need to be accounted for.
      transAt <- sum(indEqual & !halt)
      cdfAfter[i] <- (transBefore + transAt) / numSamples
    }
  }

  if (where == "before") {
    return(list(y = y.calc, before = cdfBefore))
  }

  if (where == "after") {
    return(list(y = y.calc, after = cdfAfter))
  }

  return(list(y = y.calc, before = cdfBefore, after = cdfAfter))
}
