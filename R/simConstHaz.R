#' @title Simulate event history data with constant event and censoring hazards
#'  
#' @description \code{simConstHaz} generates simulated event history data assuming the event hazard, a1, and censoring hazard, a0, are constants
#' 
# @details
#'
#' @param n Number of observations to simulate
#' @param a1 Event hazard
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
#' @seealso simLinHaz, fitConstHaz, fitLinHaz
# @references 


simConstHaz <- function(n,a1,a0) {
    eventTime     <- rexp(n,rate=a1)
    censoringTime <- rexp(n,rate=a0)
    indCensored <- eventTime > censoringTime
    transTime <- eventTime
    transTime[indCensored] <- censoringTime[indCensored]
    return(survival::Surv(transTime,!indCensored))
}

