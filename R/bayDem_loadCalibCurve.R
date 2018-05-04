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
#   calibFile    string        The calibration file
#
# Output(s)
#   Name         Type          Description
#   calibDf      data frame    T

bayDem_loadCalibCurve <- function(calibFile) {
    calibDf <- read.csv(calibFile,comment.char='#',header=F)
    calibDf <- calibDf[,1:3]
    colnames(calibDf) <- c('yearBP','uncalYearBP','uncalYearBPError')
    return(calibDf)
}
