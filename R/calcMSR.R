#' @title Calculate Martingale Survival Residuals
#'
#' @description Calculate Martingale Survival Residuals for a vector of cumulative hazards and status indicators (0/FALSE for censoring and 1/TRUE for events)
#'
# @details
#'
#' @param cumHaz A vector of cumulative hazards for the event of interest (length N)
#' @param status A vector of status indicators (length N)
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return A vector of Martingale Survival Residuals (length N)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcMSR
#' @references
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}

calcMSR <- function(cumHaz, status = NA) {
  if (all(is.na(status))) {
    # Assume all observations are events, likely because censoring has been modeled and incorporated into the cumulative hazard
    status <- rep(1, length(cumHaz))
    msr <- (-.5) * (2 * exp(-cumHaz) - 1)
  } else {
    msr <- (-.5) * ((1 + status) * exp(-cumHaz) - 1)
  }

  return(msr)
}
