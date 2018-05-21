#' @title Calculate the Negative Martingale Residual probability density
#'  
#' @description Calculate the expected probability density function for the Negative Martingale Residual. This assumes constant transition intensities for both event occurrence and censoring.
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
#' @seealso calcMartResid
#'
#' @references 
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}

calcNMRDensity <- function(y,rho,component="both") {
    f.e <- exp(-(y+1)/(1-rho)) # event component
    f.c <- rep(0,length(y))    # censoring component
    ind <- y >= 0 # Indices of censored component
    f.c[ind] <- f.c[ind] + exp(-y[ind]/(1-rho)) * rho/(1-rho)
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
