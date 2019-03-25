#' @title Do calculations for input KidStats data
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
processKidStatsInput <- function(kidStatsInput,cumProbitModel,xknown=NA,verbose=F) {
  # Extract the input data
  y <- rep(NA,length(cumProbitModel$varNames))
  if(verbose) {
    print('--')
    print('Input variable value list, kidStatsInput$varVal')
    print(kidStatsInput$varVal)
  }
  for(v in 1:length(cumProbitModel$varNames)) {
    if(cumProbitModel$varNames[v] %in% names(kidStatsInput$varVal)) {
      if(v <= cumProbitModel$hp$J) {
        # Ordinal
        originalCat = as.numeric(kidStatsInput$varVal[cumProbitModel$varNames[v]])
        yadaCat <- which(cumProbitModel$hp$originalCat[[v]] == originalCat) - 1
        y[v] <- yadaCat
        if(verbose) {
          print('-- ordinal --')
          print(cumProbitModel$varNames[v])
          print(paste('Input value:',as.character(originalCat)))
          print(paste('Remapped value for yada:',as.character(yadaCat)))
	}
      } else {
        # Continuous
          print('-- continuous --')
          print(cumProbitModel$varNames[v])
          print(as.numeric(kidStatsInput$varVal[cumProbitModel$varNames[v]]))
        y[v] <- as.numeric(kidStatsInput$varVal[cumProbitModel$varNames[v]])
      }
    }
  }

  fv <- calc_x_posterior(kidStatsInput$xcalc,y,cumProbitModel$theta_x,cumProbitModel$theta_y_list,cumProbitModel$hp)
  if(is.na(xknown)) {
    post <- analyze_x_posterior(kidStatsInput$xcalc,fv)
  } else {
    post <- analyze_x_posterior(kidStatsInput$xcalc,fv,xknown)
  }
  return(post)
}
