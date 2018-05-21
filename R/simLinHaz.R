#' @title Simulate event history data with linear event and constant censoring hazards
#'  
#' @description \code{simLinHaz} generates simulated event history data assuming the event hazard is linear, hazEvent = a1*t + b1, and the censoring hazard is constant, a0.
#' 
# @details
#'
#' @param n Number of observations to simulate
#' @param a1 Intercept for the linear event hazard
#' @param b1 Slope for the linear event hazard
#' @param a0 Censoring hazard
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


simLinHaz <- function(n,a1,b1,a0) {
    cdfEvent <- runif(n) # Take a draw of the cumulative distribution function (CDF) from 0 to 1
    eventTime    <- (a1/b1)*(sqrt(1 - 2*b1*log(1-cdfEvent)/a1^2)-1)
    censoringTime <- rexp(n,rate=a0)
    indCensored <- eventTime > censoringTime
    transTime <- eventTime
    transTime[indCensored] <- censoringTime[indCensored]
    return(data.frame(time=transTime,status=as.numeric(!indCensored)))
}

