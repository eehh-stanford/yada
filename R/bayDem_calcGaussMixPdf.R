# @keywords
#' @export
# Description
#   Calculate the probability density function (pdf) for a Gaussian mixture with
#   K components:
#
#   f = [p(y_1|th)]
#       [p(y_2|th)]
#       [   ...   ]
#       [   ...   ]
#
# Example calls(s)
#
#   f <- bayDem_calcGaussMixPdf(th,y)
# 
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vectore-like object of length 2 + 3*K with the
#                          following entries:
#                          pik   -- [K entries] Weight of the k-th mixture
#                          muk   -- [K entries] Mean of the k-th mixture
#                          sigk  -- [K entries] Standard deviation of the k-th
#                                               mixture
#   y       vector         Vector at which to calculate pdf (calendar dates)
#   doNorm  logical        A flag indicating whether to normalize to the
#                          interval of y
#
# Output(s)
#   Name    Type           Description
#   f       column vector  Output pdf
bayDem_calcGaussMixPdf <- function(th,y,doNorm=F) {
    # Since pnorm does not support fully vectorized operations of this type,
    # the inputs to pnorm must be replicated into vectors and then turned into
    # matrices
    K <- length(th)/3 # Number of  mixtures
    G <- length(y)
    yRep <- rep(y,K)
    piRep <- as.vector(t(matrix(th[1:K],G,nrow=K)))
    muRep  <- as.vector(t(matrix(th[(K+1):(2*K)],G,nrow=K)))
    sigRep <- as.vector(t(matrix(th[(2*K+1):(3*K)],G,nrow=K)))
    pdfVect <- dnorm(yRep,muRep,sigRep) * piRep

    pdfMatrix <- matrix(pdfVect,ncol=K)

    f <- rowSums(pdfMatrix)
    if(doNorm) {
      f <- f / sum(f) / (y[2]-y[1]) # assumes y is evenly spaced
    }
    return(f)
}
