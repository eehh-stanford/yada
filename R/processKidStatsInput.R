#' @title Do calculations for input KidStats data
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
processKidStatsInput <- function(kidStatsInput,cumProbitModel) {
  # Extract the input data
  y <- rep(NA,length(cumProbitModel$varNames))
  for(v in 1:length(cumProbitModel$varNames)) {
    if(cumProbitModel$varNames[v] %in% names(kidStatsInputs$varVal)) {
      y[v] <- as.numeric(kidStatsInput$varVal[cumProbitModel$varNames[v]])
    }
  }

  fv <- calc_x_posterior(kidStatsInput$xcalc,y,cumProbitModel$theta_x,cumProbitModel$theta_y_list,hp)
  post <- analyze_x_posterior(kidStatsInput$xcalc,fv,inputData$x[n])
  return(post)
}
