#' @title Calculate empirical cumulative density function point pairs for plotting
#'  
#' @description Create point pairs for plotting the input empirical cumulative density function (CDF; see \code{calcStoppingCDF} for the form of the empirical CDF)
#' 
# @details
#'
#' @param cdf An empirical CDF of the type output by calcStoppingCDF
#' 
# @keywords
#' @export
#' 
# @examples
#' 
#' @return A matrix of point pairs with each row having (x1,y1,x2,y2) for a given point pair
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>
#' @seealso calcStoppingCDF, addECDFLinesToPlot
#
#' @references 
#' \url{https://www.biorxiv.org/content/early/2017/01/31/104406}


calcECDFPointPairs <- function(cdf) {
	numRows <- 2*length(cdf$y)-1
	if(!(-.5 %in% cdf$y)) {
		numRows <- numRows + 1
		row0 <- 1 
	} else {
		row0 <- 0
	}

	if(! (.5 %in% cdf$y)) {
		numRows <- numRows + 1
	}

	# Matrix in which to store the point pair information (x1,y1,x2,y2)
	pn <- matrix(NA,nrow=numRows,ncol=4)

	if(!(-.5 %in% cdf$y)) {
		pn[1,] <- c(-.5,0,cdf$y[1],0)
	}

	N <- length(cdf$y)
	for(k in 1:N) {
		pn[(k-1)*2+row0+1,] <- c(cdf$y[k],cdf$before[k],cdf$y[k]  ,cdf$after[k])
		if(k < N) {
			pn[(k-1)*2+row0+2,] <- c(cdf$y[k],cdf$after[k] ,cdf$y[k+1],cdf$before[k+1])
		}
	}

	if(! (.5 %in% cdf$y)) {
		pn[numRows,] <- c(cdf$y[N],cdf$after[N],.5,1)
	}
	
	return(pn)
}


