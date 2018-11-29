# @keywords
#' @export
# Description
#   Simulate the radiocarbon measurement process given a set of input
#   calendar dates.
#
# Example calls(s)
#
#   sampRcMeas <- bayDem_dateSampToC14Samp(calibFile)
#
# Input(s)
#   Name          Type      Description
#   y_e           vector    The input calendar dates (_e for events)
#   calibDf       dframe    The radiocarbon calibration curve
#   errorSpec     list      Specification for measurement errors
#   isAd          boolean   Boolean indicating whether input dates are AD
#                           (default is False -- i.e., BP) [optional]
#
# Output(s)
#   Name          Type       Description
#   sampRcMeas    list       A list with the following vectors (_m for
#                            measurement):
#                            phi_m     -- C14 / C12 ratio
#                            sig_m     -- Uncertainty on phi_m
#                            yrc_m     -- Radiocarbon date from phi_m
#                            sig_yrc_m -- Uncertainty in radiocarbon years

bayDem_dateSampToC14Samp <- function(y_e, calibDf, errorSpec, isAD = F) {
  # If input calendar dates are AD, convert to calBP
  if (isAD) {
    y_e <- 1950 - y_e
  }

  N <- length(y_e)

  y_curve <- rev(calibDf$yearBP)
  mu_k_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_k_curve <- mu_k_curve * rev(calibDf$uncalYearBPError) / 8033

  # Interpolate curves at y_e to yield mu_k
  mu_k <- approx(y_curve, mu_k_curve, y_e)
  mu_k <- mu_k$y

  # Interpolate curves at y_e to yield sig_k
  sig_k <- approx(y_curve, sig_k_curve, y_e)
  sig_k <- sig_k$y

  # Sample the measurement errors
  sig_m <- runif(N, errorSpec$min, errorSpec$max)
  sig_tot <- sqrt(sig_m^2 + sig_k^2)

  # The measured "ratios"
  phi_m <- rnorm(N, mu_k, sig_tot)

  # Calculate radiocarbon years (uncal) measurement and error
  sig_yrc_m <- 8033 * sig_m / phi_m
  yrc_m <- -8033 * log(phi_m)
  return(list(phi_m = phi_m, sig_m = sig_m, yrc_m = yrc_m, sig_yrc_m = sig_yrc_m))
}
