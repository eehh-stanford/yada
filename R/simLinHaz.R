#' @title Simulate event history data with linear event and constant censoring hazards
#'
#' @description \code{simLinHaz} generates simulated event history data assuming the event hazard is linear, hazEvent = b1 + a1*t, and the censoring hazard is constant, hazCens = b2.
#'
# @details
#'
#' @param n Number of observations to simulate
#' @param b1 Intercept for the linear event hazard
#' @param a1 Slope for the linear event hazard
#' @param b2 Censoring hazard
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return A Surv object (survival package) with:
#' \describe{
#'   \item{\code{time}}{The event or censoring time}
#'   \item{\code{status}}{The status indicator: TRUE for an event and FALSE for censoring}
#' }
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz, fitConstHaz, fitLinHaz
# @references


simLinHaz <- function(n, b1, a1, b2) {
  cdfEvent <- runif(n) # Take a draw of the cumulative distribution function (CDF) from 0 to 1
  eventTime <- (b1 / a1) * (sqrt(1 - 2 * a1 * log(1 - cdfEvent) / b1^2) - 1)
  censoringTime <- rexp(n, rate = b2)
  indCensored <- eventTime > censoringTime
  transTime <- eventTime
  transTime[indCensored] <- censoringTime[indCensored]
  return(survival::Surv(transTime, !indCensored))
}
