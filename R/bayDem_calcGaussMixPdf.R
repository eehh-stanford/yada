# Description
#   Calculate the probability density function (pdf) for a (truncated)
#   Gaussian mixture with J components:
#
#   f = [p(y_1|th)]
#       [p(y_2|th)]
#       [   ...   ]
#       [   ...   ]
#
# Example calls(s)
#
#   bayDem_calcGaussMixPdf <- function(th,y) {
# 
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vectore-like object of length 2 + 3*J with the
#                          following entries:
#                          ymin  -- Minimum value
#                          ymax  -- Maximum value
#                                   [Samples are truncated to the interval
#                                    ymin to ymax]
#                          zj    -- [J entries] Weight of the j-th mixture
#                          muj   -- [J entries] Mean of the j-th mixture
#                          sigj  -- [J entries] Standard deviation of the j-th
#                                               mixture
#   y       vector         Vector at which to calculate pdf (calendar dates)
#
# Output(s)
#   Name    Type           Description
#   f       column vector  Output pdf
bayDem_calcGaussMixPdf <- function(th,y) {
    J <- (length(th)-2)/3 # Number of  mixtures
    # Pre-calculate the normalization vector that accounts for truncation
    normFact <- sum((pnorm(th[2],th[(3+J):(2+2*J)],th[(3+2*J):(2+3*J)]) - pnorm(th[1],th[(3+J):(2+2*J)],th[(3+2*J):(2+3*J)]))*th[3:(2+J)])

    # Combine the mixture proportions with the normalization factor for the overall weighting
    weightVect <- th[3:(2+J)] / normFact
    
    f <- matrix(0,length(y),1)


    # Create a matrix of the pdf values of size [N' x J], where N' is the
    # number entries in y that are inside ymin and ymax. Since pnorm does not
    # support fully vectorized operations of this type, the inputs to pnorm
    # must be replicated into vectors and then turned into matrices
    # B is true for non-truncated samples
    B <- th[1] <= y & y <= th[2]
    yRep <- rep(y[B],J)
    muRep <- as.vector(t(matrix(th[(3+J):(2+2*J)],sum(B),nrow=J)))
    sigRep <- as.vector(t(matrix(th[(3+2*J):(2+3*J)],sum(B),nrow=J)))
    pdfVect <- dnorm(yRep,muRep,sigRep)

    pdfMatrix <- matrix(pdfVect,ncol=J) # N' x J

    # Calculate the mixture PDF, accounting for the overall weighting of each
    # component
    weightMatrix <- t(matrix(rep(weightVect,sum(B)),nrow=J))

    f <- rep(0,length(y))
    f[B] <- rowSums(pdfMatrix * weightMatrix)
    return(f)
}
