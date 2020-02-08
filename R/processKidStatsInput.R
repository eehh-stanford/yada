#' @title Do calculations for input KidStats data
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
processKidStatsInput <- function(kidStatsInput,cumProbitModel,catList,xknown=NA) {
  # Extract the input data
  y <- y_from_ks(kidStatsInput,cumProbitModel,catList)

  fv <- calc_x_posterior(kidStatsInput$xcalc,y,cumProbitModel$theta_x,cumProbitModel$theta_y_list,cumProbitModel$hp)
  if(is.na(xknown)) {
    post <- analyze_x_posterior(kidStatsInput$xcalc,fv)
  } else {
    post <- analyze_x_posterior(kidStatsInput$xcalc,fv,xknown)
  }
  return(post)
}

y_from_ks <- function(kidStatsInput,model,catList) {
  y <- rep(NA,length(model$varNames))

  for(vv in 1:length(kidStatsInput$varVal)) { 
    varName <- names(kidStatsInput$varVal)[vv]
    varIndex <- which(model$varNames == varName)
    if(varIndex <= model$hp$J) {
      # A boolean vector
      B <- catList[[varName]] == kidStatsInput$varVal[varName]
      if(sum(B) != 1) {
        stop(paste('Wrong number of matches for variable',varName))
      }
      y[varIndex] <- which(B) - 1
    } else {
      y[varIndex] <- kidStatsInput$varVal[[vv]]
    }
  }
  return(y)
}
