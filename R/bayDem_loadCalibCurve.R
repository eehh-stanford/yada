#' @title Load Calibration Curve
#'  
#' @description Parse and return the radiocarbon calibration curve stored in data
#' 
# @details
#'
#' @param calibCurve Name of calibration curve
#' @param status A vector of status indicators (length N)
#' 
# @keywords
#' @export
#' 
# @examples
#' 
#' @return The calibration dataframe, with columns yearBP, uncalYearBP, and uncalYearBPError
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

bayDem_loadCalibCurve <- function(calibCurve) {
  if(calibCurve == 'intcal13') {
    calibFile <- system.file('data/intcal13.14c.data',package='yada')
  } else {
    stop(paste('Unknown calibration curve name:',calibCurve))
  }
    calibDf <- read.csv(calibFile,comment.char='#',header=F)
    calibDf <- calibDf[,1:3]
    colnames(calibDf) <- c('yearBP','uncalYearBP','uncalYearBPError')
    return(calibDf)
}
