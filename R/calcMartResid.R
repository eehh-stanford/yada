#' @title Calculate Martingale Residuals
#'
#' @description Calculate Martingale Residuals for a vector of cumulative hazards and status indicators (0/FALSE for censoring and 1/TRUE for events)
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
#' @return A vector of Martingale Residuals (length N)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcMSR
#' @references
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}

calcMartResid <- function(cumHaz, status) {
  return(status - cumHaz)
}
