#' @title Calculate the Martingale Residual probability density
#'  
#' @description Calculate the expected probability density function for the Martingale Residual. This assumes constant transition intensities for both event occurrence and censoring.
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

calcMartResidDensity <- function(y,rho,component="both") {
    return(calcNMRDensity(-y,rho,component))
}
