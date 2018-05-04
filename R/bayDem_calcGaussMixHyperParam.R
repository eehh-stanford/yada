# Description
#   Calculate and return the (truncated) Gaussain mixture hyper parameter
#   object.
#
#   (1) The mean (mu) is drawn uniformly from ymin to ymax
#
#   (2) The scale (sigma) is drawn from an offset gamma distribution:
#
#       sigma ~ sig_startOffset + gamma(shape=sig_alpha,rate=sig_beta)
#
#       For a gamma distribution, the following holds for the mode [sig_mode0
#       b/c momentarily ignoring offset]:
#
#       (a) sig_mode0 = (sig_alpha-1)/beta
#           --> sig_beta = (sig_alpha-1)/sig_mode0
#
#       The input location of the mode, sig_modeOffset, is relative to the
#       origin, so
#
#       sig_mode0 = sig_modeOffset - sig_startOffset
#
#       Hence:
#
#       sig_beta = (sig_alpha-1)/(sig_modeOffset - sig_startOffset)
#
#   (3) The mixture proportions are drawn from the dirichlet distribution with
#       parameter z_alpha -> rep(z_alpha,J), where J is the number of mixtures.
#
# Example calls(s)
#
#   gmPriorParam <- bayDem_calcGaussMixPriorParam(ymin,ymax,startOffset,modeOffset,gammaAlpha,dirichParam)
# 
# Input(s)
#   Name            Type      Description
#   ymin            scalar    Minimum date
#   ymax            scalar    Maximum date
#   sig_startOffset scalar    Constant added to gamma draw for the scale (sigma)
#   sig_modeOffset  scalar    Location of the mode for gamma draw for the scale
#                             [see above]
#   sig_alpha       scalar    alpha for the gamma draw for the sigma
#   z_alpha         vector    dirichlet parameter vector for the mixture
#                             proportions
#
# Output(s)
#   Name       Type     Description
#   gmHyper    list     The hyperparameters:
#                         $ymin          -- Minimum calendar date
#                         $ymax          -- Maximum calendar date
#                         $z$zAlpha      -- Parameter for the Dirichlet draw
#                         $sig$sigAlpha  -- Shape parameter for the gamma draw
#                         $sig$sigBeta   -- Rate parameter for the gamma draw
#                         $sig$sigOffset -- Offset for the gamma draw
#

bayDem_calcGaussMixHyperParam <- function(ymin,ymax,sig_startOffset,sig_modeOffset,sig_alpha,z_alpha) {
    gmHyper <- list(ymin=ymin,ymax=ymax)
    gmHyper$z$zAlpha <- z_alpha

    sig_beta  <- (sig_alpha-1)/(sig_modeOffset-sig_startOffset)
    gmHyper$sig$sigAlpha <- sig_alpha
    gmHyper$sig$sigBeta <- sig_beta
    return(gmHyper)
}

