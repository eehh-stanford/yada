#' @title Combinadic functions
#'
#' @description Implementation of combinadic-related functions as described here:
#'
#' James McCaffrey -- Generating the mth Lexicographical Element of a Mathematical Combination
#'
#' Indexing begins from 0, not 1.
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export
comb <- function(m,N,k) {
  # The combinadic
  x <- m
  a <- N
  b <- k

  # Use cv instead of c since c is an R function
  cv = rep(NA,k)

  for(ii in 1:k) {
    cv[ii] <- largestV(a,b,x)
    x <- x - nchoosek_smart(cv[ii],b)
    a <- cv[ii]
    b <- b - 1
  }
  return(cv)
}

#' @export
dual <- function(m,nElements=NA,N=NA,k=NA) {
  # The dual
  if(!is.na(nElements)) {
    d <- (nElements-1) - m
    return(d)
  }

  if(!is.na(N) && !is.na(k)) {
    d <- (nchoosek_smart(N,k)-1) - m
    return(d)
  }

  stop('Unrecognized input pattern')
}

#' @export
elem <- function(m,N,k) {
  # The element
  x <- dual(m,N=N,k=k)
  
  # Use cv instead of c since c is an R function
  cv <- comb(x,N,k)
  e <- (N-1) - cv
  return(e)
}

#' @export
elemToIndex <- function(e,N) {
  # Convert the element to an index
  k <- length(e)
  y <- 0
  for(ii in 1:k) {
    y <- y + nchoosek_smart(N-1-e[ii],k-(ii-1));
  }

  index = nchoosek_smart(N,k) - 1 - y
  return(index)
}

#' @export
largestV <- function(a,b,x) {
  # Calculate the largest V
  v <- a-1
  while(nchoosek_smart(v,b) > x) {
    v <- v - 1;
  }
  return(v)
}

#' @export
nchoosek_smart <- function(N,k) {
  # A smart calculation of n choose k

  if(N < k) {
    nck <- 0
  } else if(N==k) {
    nck <- 1
  } else {
    nck <- choose(N,k)
  }
  return(nck)
}
