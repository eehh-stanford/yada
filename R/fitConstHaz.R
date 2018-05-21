#' @title Fit event history data assuming constant event and censoring hazards
#'  
#' @description \code{fitConstHaz} does a maximum likelihood fit of the input data assuming constant event and censoring hazards. The input is a Surv object (package survival) with a time column and a status column. Time is the event or censoring time. Status is 1 (or TRUE) if the event occured and 0 (or FALSE) if censoring occurred.
#' 
# @details
#'
#' @param survObj A Surv object (survival package) with event times and a status indicator (see \code{simConstHaz})
#' 
# @keywords
#' @export
#' 
# @examples
#' 
#' @return A list with three fields
#'
#' \describe{
#'   \item{\code{a}}{Overall, combined event and censoring rate}
#'   \item{\code{a1}}{Event rate}
#'   \item{\code{a0}}{Censoring rate}
#' }
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso simConstHaz
# @references 

fitConstHaz <- function(survObj) {
    a <- 1/mean(survObj[,'time'])
    a1 <- a * sum(survObj[,'status']) / length(survObj[,'status'])
    a0 <- a * sum(!survObj[,'status']) / length(survObj[,'status'])
    return(list(a=a,a1=a1,a0=a0))
}
