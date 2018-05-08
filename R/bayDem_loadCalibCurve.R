# Description
#   Parse and return the radiocarbon calibration curve stored in the input file
#   calibFile.
#
# Example calls(s)
#
#   calibDf <- bayDem_loadCalibCurve(calibFile)
# 
# Input(s)
#   Name         Type          Description
#   calibCurve   string        Name of calibration curve
#
# Output(s)
#   Name         Type          Description
#   calibDf      dataframe     The calibration dataframe, with columns yearBP,
#                              uncalYearBP, and uncalYearBPError

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
