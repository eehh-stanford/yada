#' @title Calculate the Martingale Survival Residual probability density
#'  
#' @description Calculate the expected probability density function for the Martingale Survival Residual (MSR). This assumes constant transition intensities for both event occurrence and censoring.
#' 
# @details
#'
#' @param y A vector of locations at which to calculate the probability density
#' @param rho The censoring ratio
#' @param component The component of the probability density to return. "both" yields the full density. "event" yields the component due to event occurrence. "censored" yields the component due to censoring. The default is "both".
#' 
# @keywords
#' @export
#' 
# @examples
#' 
#' @return A vector of probability densities that is the same length as y
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcMSR
#'
#' @references 
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}

calcMSRDensity <- function(y,rho,component="both") {
    f.e <- (.5-y)^(rho/(1-rho)) # event component
    f.c <- rep(0,length(y)) # censoring component
    ind <- y >= 0 # Indices of censored component
    f.c[ind] <-  f.c[ind] + (1-2*y[ind])^(rho/(1-rho))*2*rho/(1-rho)
    f <- f.e + f.c
    if(component=="event") {
        return(f.e)
    } else if(component=="censored") {
        return(f.c)
    } else if(component=="both") {
        return(f)
    } else {
        stop(paste("Unrecognized component argument:",component))
    }
}
