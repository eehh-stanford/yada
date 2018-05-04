# Description
#   Sample from a (truncated) two-component Gaussian mixture. This provides
#   calendar dates of radiocarbon samples from the demographic model specified
#   by the two-component Gaussian mixture.
#
# Example calls(s)
#
#   samp <- bayDem_sampleGaussMix(N,th)
# 
# Input(s)
#   Name    Type           Description
#   N       integer        The number of samples
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
#
# Output(s)
#   Name    Type           Description
#   samp    vector         The samples (length = N)

bayDem_sampleGaussMix <- function(N,th) {
    J <- (length(th)-2)/3 # Number of  mixtures

    # This commented out code fails because Truncate expects an
    # AbscontDistribution but creating a mixture distribution via vectors
    # as in this commented out code produces a UnivarMixingDistribution, which
    # cannot be cast to an AbscontDistribution. Perhaps this shortcoming of
    # distr will be addressed in the future. In the meantime, the number of
    # mixtures is limited to 4

    #mixDistrList <- list()
    #for(j in 1:J) {
    #    mixDistrList[[j]] <- distr::Norm(mean=th[2+J+j],sd=th[2+2*J+j])
    #}
    #mixDistr <- new('UnivarMixingDistribution',mixCoeff=th[3:(2+J)],mixDistr=new('UnivarDistrList',mixDistrList))
    #mixDistrTrunc <- distr::Truncate(mixDistr,th[1],th[2])
    
    # (Somewhat awkwardly), define the mixing distribution directly for J = 1
    # to 6 directly
    if(J == 1) { # Allow J = 1 (i.e., a Gaussian, not a mixture) as a special
	         # case
        normMix <- distr::Norm(mean=th[4],sd=th[5])
    } else if(J == 2) {
        normMix <- distr::UnivarMixingDistribution(distr::Norm(mean=th[5],sd=th[7]),distr::Norm(mean=th[6],sd=th[8]),mixCoeff=th[3:4])
    } else if(J == 3) {
        normMix <- distr::UnivarMixingDistribution(distr::Norm(mean=th[6],sd=th[9]),distr::Norm(mean=th[7],sd=th[10]),,distr::Norm(mean=th[8],sd=th[11]),mixCoeff=th[3:5])
    } else if(J == 4) {
        normMix <- distr::UnivarMixingDistribution(distr::Norm(mean=th[7],sd=th[11]),distr::Norm(mean=th[8],sd=th[12]),,distr::Norm(mean=th[9],sd=th[13]),distr::Norm(mean=th[10],sd=th[14]),mixCoeff=th[3:6])
    } else {
        stop('The maximum number of supported mixture components is 4')
    }

    # Use distr to sample from a two-component, truncated Gaussian mixture
    normMixTrunc <- distr::Truncate(normMix,th[1],th[2])
    samp <- distr::r(normMixTrunc)(N)
    return(samp)
}

