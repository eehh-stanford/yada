# Description
#   Calculate the probability density function (pdf) for a (truncated)
#   Gaussian mixture with K components:
#
#   f = [p(y_1|th)]
#       [p(y_2|th)]
#       [   ...   ]
#       [   ...   ]
#
# Example calls(s)
#
#   bayDem_calcGaussMixPdf <- function(th,y)
# 
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vectore-like object of length 2 + 3*K with the
#                          following entries:
#                          ymin  -- Minimum value
#                          ymax  -- Maximum value
#                                   [Samples are truncated to the interval
#                                    ymin to ymax]
#                          pik   -- [K entries] Weight of the k-th mixture
#                          muk   -- [K entries] Mean of the k-th mixture
#                          sigk  -- [K entries] Standard deviation of the k-th
#                                               mixture
#   y       vector         Vector at which to calculate pdf (calendar dates)
#
# Output(s)
#   Name    Type           Description
#   f       column vector  Output pdf
bayDem_calcGaussMixPdf <- function(th,y) {
    K <- (length(th)-2)/3 # Number of  mixtures
    # Pre-calculate the normalization vector that accounts for truncation
    normFact <- sum((pnorm(th[2],th[(3+K):(2+2*K)],th[(3+2*K):(2+3*K)]) - pnorm(th[1],th[(3+K):(2+2*K)],th[(3+2*K):(2+3*K)]))*th[3:(2+K)])

    # Combine the mixture proportions with the normalization factor for the overall weighting
    weightVect <- th[3:(2+K)] / normFact
    
    f <- matrix(0,length(y),1)


    # Create a matrix of the pdf values of size [N' x K], where N' is the
    # number entries in y that are inside ymin and ymax. Since pnorm does not
    # support fully vectorized operations of this type, the inputs to pnorm
    # must be replicated into vectors and then turned into matrices
    # B is true for non-truncated samples
    B <- th[1] <= y & y <= th[2]
    yRep <- rep(y[B],K)
    muRep <- as.vector(t(matrix(th[(3+K):(2+2*K)],sum(B),nrow=K)))
    sigRep <- as.vector(t(matrix(th[(3+2*K):(2+3*K)],sum(B),nrow=K)))
    pdfVect <- dnorm(yRep,muRep,sigRep)

    pdfMatrix <- matrix(pdfVect,ncol=K) # N' x K

    # Calculate the mixture PDF, accounting for the overall weighting of each
    # component
    weightMatrix <- t(matrix(rep(weightVect,sum(B)),nrow=K))

    f <- rep(0,length(y))
    f[B] <- rowSums(pdfMatrix * weightMatrix)
    return(f)
}
