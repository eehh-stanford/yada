# Description
#   Calculate the probability density for the signal parametrized by augData,
#   the augmented data. augData accounts for coefficients that are zero via the
#   indicator matrix S.
#
# Example calls(s)
#
#   f <- bayDem_calcWaveletPdf(augData)
#   f <- bayDem_calcWaveletPdf(augData,5)
# 
# Input(s)
#   Name       Type      Description
#   augData    list      The augmented data
#
# Output(s)
#   Name       Type      Description
#   f          vector    The probability density function 

bayDem_calcWaveletPdf <- function(augData,dy=1) {
    J <- length(augData$S)
    dwt <- augData$dwt
    # Iterate over S to place zeros in augData$m as appropriate
    for(j in 1:J) {
        dwt$data[[j]] <- dwt$data[[j]] * augData$S[[j]]
    }
    f_root <- reconstruct(dwt)
    f <- f_root^2
    f <- f / sum(f) / dy
    return(f)
}
