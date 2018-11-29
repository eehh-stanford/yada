#' @title Do Demographic Bayesian Inference
#'
#' @description This is the core function that implements the Bayesian inference. The input is a problem statement object (a list), prob, that consists of the input data (the vectors phi_m and sig_m) and the hyperparameters (hp). stan is called via the rstan package to sample from the posterior. The output is the variable soln of class bayDem_soln, which is a list with the fields prob (the input) and fit (the result of the stan fit). prob can also have an optional field control that specifies the following control parameters for the Bayesian inference [default in brackets]:
#'   numChains     -- [4]    Number of chains
#'   sampsPerChain -- [2000] Number of samples per chain
#'   initList      --        The initializations for each chain. The default is
#'                           to sample from prior using hyperparameters
#'
# @details
#' @param prob List with the fields phi_m (vector of radiocarbon measurements as fraction modern), sig_m (vector of measurement errors for phi_m), hp (list of hyperparameters), and calibDf (calibration dataframe). In addition, the field control is optional (see above).
#'
# @keywords
#' @export
#'
# @examples
#'
#' @return soln, a list with three fields: prob (the input), fit (the result of the call to stan), and control (the control parameters used)
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

bayDem_doInference <- function(prob) {
  # Unpack and/or define the control parameters
  if (exists("control", where = prob) == T) {
    haveNumChains <- exists("numChains", where = prob$control) == T
    haveSampsPerChain <- exists("sampsPerChain", where = prob$control) == T
    haveInitList <- exists("initList", where = prob$control) == T
    haveStanControl <- exists("stanControl", where = prob$control) == T
  } else {
    haveNumChains <- F
    haveSampsPerChain <- F
    haveInitList <- F
  }

  if (haveNumChains) {
    numChains <- prob$control$numChains
  } else {
    numChains <- 4
  }

  if (haveSampsPerChain) {
    sampsPerChain <- prob$control$sampsPerChain
  } else {
    sampsPerChain <- 2000
  }

  if (haveInitList) {
    initList <- prob$control$initList
  } else {
    initList <- bayDem_samplePrior(prob$hp, numChains)
  }

  if (haveStanControl) {
    stanControl <- prob$control$stanControl
  } else {
    stanControl <- NA
  }

  controlFinal <- list(numChains = numChains, sampsPerChain = sampsPerChain, initList = initList, stanControl = stanControl)

  if (prob$hp$fitType == "gaussmix") {
    # Stan needs all the inputs and hyperparameters as variables in R's workspace
    ymin <- prob$hp$ymin
    ymax <- prob$hp$ymax
    ygrid <- seq(ymin, ymax, by = prob$hp$dy)
    M <- bayDem_calcMeasMatrix(ygrid, prob$phi_m, prob$sig_m, calibDf)
    Mt <- t(M)
    N <- dim(M)[1]
    G <- dim(M)[2]
    sigAlpha <- prob$hp$sigAlpha
    sigBeta <- prob$hp$sigBeta
    dirichParam <- prob$hp$dirichParam
    K <- prob$hp$K
    filePath <- system.file("stan/gaussmix.stan",
      package = "yada"
    )
    options(mc.cores = parallel::detectCores())
    if (haveStanControl) {
      fit <- stan(filePath, chains = numChains, iter = sampsPerChain, init = initList, control = stanControl)
    } else {
      fit <- stan(filePath, chains = numChains, iter = sampsPerChain, init = initList)
    }
  } else {
    stop(paste("Unrecognized fit type:", prob$hp$fitType))
  }

  soln <- list(prob = prob, fit = fit, control = controlFinal)
  class(soln) <- "bayDem_soln"
  return(soln)
}
