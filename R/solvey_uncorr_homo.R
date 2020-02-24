#' @title Do a maximum likelihood fit for theta_y for the uncorrelated, homoskedastic case
#'
#' @description x is a vector of length N, Y a matrix of dimensions (J+K) by N,
#' and hp a list of hyperparameters. For each variable j (ordinal) and k
#' (continuous), do a maximum likelihood fit. These fits are entirely
#' independent for the uncorrelated, homoskedastic case.
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
solvey_uncorr_homo <- function(x,Y,hp,verbose=F) {

  # The number of ordinal variables
  J <- hp$J
  # The number of continuous variables
  K <- hp$K

  # Initialize the fit

  if(J > 0) {
    rho <- rep(NA,J)
    tau <- list()
  }

  if(K > 0) {
    a <- rep(NA,K)
    r <- rep(NA,K)
    b <- rep(NA,K)
  }
  s  <- rep(NA,J+K)

 if(J > 0) {

  ordFitList <- list()
  for(j in 1:J) {
    if(verbose) {
      print(paste('j =',j,'of',J))
    }
    xj <- x
    vj <- Y[j,]
    xj <- xj[!is.na(vj)]
    vj <- vj[!is.na(vj)]
    hp_j <- list(paramModel='powLawOrdHomo')
    hp_j$J <- 1
    hp_j$M <- length(unique(vj)) - 1
    ordFit <- fitPowLawOrd(xj,vj,hetero=F)
    ordFitList[[j]] <- ordFit

    th_v  <- ordFit$th_v
    th_v0 <- ordFit$th_v0
    th_v_list <- theta_y_vect2list(th_v,hp_j)
    rho[j] <- th_v_list$rho
    tau[[j]] <- th_v_list$tau[[1]]
    s[j] <- th_v_list$s
  }
 }

 if(K > 0) {
  contFitList <- list()
  for(k in 1:K) {
    if(verbose) {
      print(paste('k =',k,'of',K))
    }
    xk <- x
    wk <- Y[J+k,]
    xk <- xk[!is.na(wk)]
    wk <- wk[!is.na(wk)]
    contFit <- fitPowLaw(xk,wk,hetero=F)
    contFitList[[k]] <- contFit
    # th_w has ordering [a,r,b,s]
    th_w <- contFit
    a[k] <- th_w[1]
    r[k] <- th_w[2]
    b[k] <- th_w[3]
    s[J+k] <- th_w[4]
  }
 }

  theta_y_list <- list(paramModel='powLawMixUncorrHomo')
  if(J > 0) {
    theta_y_list$rho <- rho
  }

  if(K > 0) {
    theta_y_list$a   <- a
    theta_y_list$r   <- r
    theta_y_list$b   <- b
  }

  if(J > 0) {
    theta_y_list$tau <- tau
  }

  theta_y_list$Sigma <- diag(s^2)

  return(theta_y_list)
}
