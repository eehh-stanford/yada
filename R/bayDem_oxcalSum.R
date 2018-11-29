# Description
#   Sum the result of a call to oxcalCalibrate
#
# Example calls(s)
#
#   M <- bayDem_oxcalSum(ygrid,phi_m,sig_m,calibDf)
#
# Input(s)
#   Name          Type      Description
#   calibObj      list      A list of class oxcAARCalibratedDatesList
#
# Output(s)
#   Name          Type      Description
#   spdf          list      A list with the summed probabilities and
#                           associated calendar dates

bayDem_oxcalSum <- function(calibObj) {
  # First, get the range of calendar dates
  y <- c()
  for (ii in 1:length(calibObj)) {
    calibi <- calibObj[[ii]]$raw_probabilities
    y <- sort(unique(c(y, calibi$dates)))
  }
  if (length(unique(diff(y))) != 1) {
    stop("Calendar dates are not a regular grid")
  }
  f_spdf <- rep(0, length(y))
  for (ii in 1:length(calibObj)) {
    calibi <- calibObj[[ii]]$raw_probabilities
    for (jj in 1:length(calibi$probabilities)) {
      ind <- which(calibi$dates[jj] == y)
      f_spdf[ind] <- f_spdf[ind] + calibi$probabilities[jj]
    }
  }
  f_spdf <- f_spdf / sum(f_spdf) / (y[2] - y[1])
  return(list(y = y, f_spdf = f_spdf))
}
