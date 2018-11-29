# @keywords
#' @export
# Description
#   Calculate the measurement matrix, which is the likelihood of the radiocarbon
#   measurements phi_m with measurement uncertainty sig_m calculated at the
#   calendar dates in the vector ygrid:
#
#  M = [p(phi_m,1|y_1)  p(phi_m,1:y_2) ... ]
#    = [p(phi_m,2|y_1)       ---       ... ]
#    = [     ...             ---       ... ]
#
#   The total uncertainty of the measurement comes from measurement error
#   (SIG_M, calculated using the measurement error for each measurement) and
#   the calibration curve error (SIK_k, calculated using the uncertainty for
#   the calibration curve at each grid point). These uncertainties (and the
#   associated measurements) should already be "projected" to 1950 equivalents
#   (e.g., in bayDem_dateSampToC14Samp).
#
# Example calls(s)
#
#   M <- bayDem_calcMeasMatrix(ygrid,phi_m,sig_m,calibDf)
#
# Input(s)
#   Name          Type      Description
#   ygrid         vector    The locations at which to calculate the likelihood
#   phi_m         vector    The radiocarbon measurements (fraction modern)
#   sig_m         vector    Uncertainty of radiocarbon measurements
#   calibDf       dframe    Calibration curve (see bayDem_loadCalibCurve)
#
# Output(s)
#   Name          Type      Description
#   M             matrix    [Nmeas x Ngrid] The measurement matrix -- that is,
#                           the likelihood of this measurement calculated for
#                           the calendar dates in ygrid.

bayDem_calcMeasMatrix <- function(ygrid, phi_m, sig_m, calibDf, normRows = T, addCalibUnc = T) {
  # ygrid is in AD
  ygrid_BP <- 1950 - ygrid

  # calibration curve
  y_curve <- rev(calibDf$yearBP)
  mu_k_curve <- exp(-rev(calibDf$uncalYearBP) / 8033)
  sig_k_curve <- rev(calibDf$uncalYearBPError) * mu_k_curve / 8033

  # Interpolate curves at ygrid_BP to yield mu_k and sig_k
  mu_k <- approx(y_curve, mu_k_curve, ygrid_BP)
  mu_k <- mu_k$y
  sig_k <- approx(y_curve, sig_k_curve, ygrid_BP)
  sig_k <- sig_k$y

  PHI_m <- replicate(length(ygrid_BP), phi_m)
  SIG_m <- replicate(length(ygrid_BP), sig_m)

  MU_k <- t(replicate(length(phi_m), mu_k))
  if (addCalibUnc) {
    SIG_k <- t(replicate(length(sig_m), sig_k))
    SIG_sq <- SIG_m^2 + SIG_k^2
  } else {
    SIG_sq <- SIG_m^2
  }

  M <- exp(-(PHI_m - MU_k)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2 * pi)

  # Add the integration widths
  G <- length(ygrid)
  dyVect <- rep(NA, length(ygrid))
  indCent <- 2:(G - 1)
  dyVect[indCent] <- (ygrid[indCent + 1] - ygrid[indCent - 1]) / 2
  dyVect[1] <- (ygrid[2] - ygrid[1]) / 2
  dyVect[G] <- (ygrid[G] - ygrid[G - 1]) / 2
  dyMat <- t(replicate(length(phi_m), dyVect))
  M <- M * dyMat

  if (normRows) {
    M <- M * matrix(1 / rowSums(M), dim(M)[1], dim(M)[2]) # Normalize the rows to sum to 1
  }
  return(M)
}
