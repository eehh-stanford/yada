#' @title Combinadic functions
#'
#' @description Implementation of combinadic-related functions as described here:
#'
#' James McCaffrey -- Generating the mth Lexicographical Element of a Mathematical Combination
#'
#' Indexing begins from 0, not 1.
#'
#' @details
#'
#' Consider a set of N things of which k unique ones are chosen. The number of
#' unique sets of things that can be chosen is choose(N,k). For example, if N=5
#' and k=3 there are choose(5,3) = 10 unique such sets. Adopting lexigraphic
#' ordering starting at 0, the combinadic (comb), dual, and element (elem) are:
#'
#' m	comb	dual	elem
#' 0	210	9	012
#' 1	310	8	013
#' 2	320	7	014
#' 3	321	6	023
#' 4	410	5	024
#' 5	420	4	034
#' 6	421	3	123
#' 7	430	2	124
#' 8	431	1	134
#' 9	432	0	234
#'
#' If the lexical index, m, is 2, the corresponding element with N=5 and k=3 is:
#'
#' e <- elem(2,5,3)
#' print(e)
#' [1] 0 1 4
#'
#' If the element is 421, the corresponding lexical index is:
#'
#' m <- elemToIndex(c(1,2,3),5)
#' print(m)
#' [1] 6
#' 
#' @param N The number of unique things
#' @param k The number of things to choose out of N total things
#' @param m The lexical index of each unique choice (starting from 0)
#' @param c The combinadic (a vector; see details)
#' @param d The dual(see details)
#' @param e The element (a vector; see details)
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
