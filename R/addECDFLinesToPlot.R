#' @title Add lines for the empirical cumulative density function to a plot
#'
#' @description Point pairs for plotting an empirical cumulative density function (CDF; see \code{calcStoppingCDF} are created by \code{calcECDFPointPairs}. Add lines representing these point pairs to a plot.
#'
# @details
#'
#' @param pn An empirical CDF of the type output by calcStoppingCDF
#' @param lwParam A parameter
#'
# @keywords
#' @export
#'
# @examples
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcStoppingCDF, calcCDFPointPairs
#
#' @references
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}


addECDFLinesToPlot <- function(pn, lwParam) {
  numRows <- dim(pn)[1]
  for (row in 1:numRows) {
    lines(pn[row, c(1, 3)], pn[row, c(2, 4)], lwd = lwParam)
  }
}
