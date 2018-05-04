# Description
#   Calculate the probability density function (pdf) for a (truncated)
#   two-component Gaussian mixture at the locations in the vector y.
#
# Example calls(s)
#
#   bayDem_calcGaussMixPdf <- function(th,y) {
# 
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vectore-like object with the following entries:
#                          (1) z1    -- Weight of the first mixture
#                          (2) z2    -- Weight of the second mixture
#                          (3) mu1   -- Mean of the first mixture
#                          (4) sig1  -- Standard deviation of the first mixture
#                          (5) mu2   -- Mean of the second mixture
#                          (6) sig2  -- Standard deviation of the second mixture
#                          (7) ymin  -- Minimum value
#                          (8) ymax  -- Maximum value
#                                       [Samples are truncated to the interval
#                                        ymin to ymax]
#   y       vector         Vector at which to calculate pdf (calendar dates)
#
# Output(s)
#   Name    Type           Description
#   pdfMat  vector         Output pdf matrix [Ny x Nmc]


bayDem_calcGaussMixPdfMat <- function(th,y) {
	Nmc <- dim(th)[2]
	Ny  <- length(y)
	
	pdfMat <- matrix(0,Ny,Nmc)

	# Pre-calculate the normalization vector that accounts for truncation
	normVect1 <- pnorm(th['ymax',],th['mu1',],th['sig1',]) - pnorm(th['ymin',],th['mu1',],th['sig1',])
	normVect2 <- pnorm(th['ymax',],th['mu2',],th['sig2',]) - pnorm(th['ymin',],th['mu2',],th['sig2',])

	for(ii in 1:Ny) {
		f1 <- rep(0,Nmc)
		f2 <- rep(0,Nmc)
		# B is true for non-truncated samples
		B <- th['ymin',] <= y[ii] & y[ii] <= th['ymax',]
		f1[B] <- dnorm(y[ii],th['mu1',B],th['sig1',B]) / normVect1[B]
		f2[B] <- dnorm(y[ii],th['mu2',B],th['sig2',B]) / normVect2[B]
		pdfMat[ii,] <- th['z1',] * f1 + th['z2',] * f2
	}
	return(pdfMat)
}
