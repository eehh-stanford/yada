#' @title Do a maximum likelihood fit for theta_y for the uncorrelated, heteroskedastic case
#'
#' @description x is a vector of length N, Y a matrix of dimensions (J+K) by N,
#' and hp a list of hyperparameters. For each variable j (ordinal) and k
#' (continuous), do a maximum likelihood fit.
#'
#' @param x A vector of independent variables [length N]
#' @param Y A matrix of ordinal and/or continuous response variables [(J+K) by N]
#' @param hp A list of hyperparameters
#' @param verbose Whether to print out information as the optimization proceeds (default FALSE)
#'
#' @return The fit as a list, theta_y_list
#'
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' @export

solvey_uncorr_hetero <- function(x,Y,hp,verbose=F,th_y0=NA) {

  hp <- problem$hp
  hp$paramModel <- 'powLawMixUncorrHetero'
  if(all(is.na(th_y0))) {
    # Initialize with the homoskedastic case
    th_y_list_homo0 <- solvey_uncorr_homo(x,Y,hp,verbose)
    th_y_list0 <- th_y_list_homo0
    th_y_list0$paramModel <- 'powLawMixUncorrHetero'
    th_y_list0$kap <- 0.0001
    th_y0 <- theta_y_list2vect(th_y_list0)
  }

  th_y_bar0 <- theta_y_constr2unconstr(th_y0,hp)

  optimControl <- list(reltol=1e-12,maxit=10000000,trace=100)
  fit <- optim(th_y_bar0,powLawMixNegLogLik,control=optimControl,x=x,Y=Y,hp=hp,transformVar=T,method='BFGS')

  th_y <- theta_y_unconstr2constr(fit$par,hp)
  th_y_list <- theta_y_vect2list(th_y,hp)
  # Both s and Sigma are included
  th_y_list$Sigma <- diag(th_y_list$s^2)
  th_y_list$paramModel <- "GenCRRA_uncorr_heterosk"
  return(th_y_list)
}
