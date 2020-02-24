# Test the functions in powLawMix

# From random.org between 1 and 1,000,000
set.seed(937974)

# Simulate data with two ordinal and two continuous variables with no
# correlations
if(F) {
kappa <- .01
rho1 <- .75
rho2 <- .50
tau1 <- c(2,5)
tau2 <- c(1,3)
s1   <- .01
s2   <- .02
th_vo1 <- c(rho1,tau1,s1)
th_vo2 <- c(rho2,tau2,s2)
hp_vo <- list(paramModel='powLawOrdUncorrHomo')
hp_vo$J = 1
hp_vo$K = 0
hp_vo$M <- 2
#powLawOrdNegLogLik(
}

#stophere
th_y_sim_homo <- list(paramModel='powLawMixUncorrHomo')
th_y_sim_homo$rho <- c(.75,.5)
th_y_sim_homo$tau <- list()
th_y_sim_homo$tau[[1]] <- c(2,5)
th_y_sim_homo$tau[[2]] <- c(1,3)
th_y_sim_homo$a <- c(2,3)
th_y_sim_homo$r <- c(.45,.10)
th_y_sim_homo$b <- c(1.2,-.5)
th_y_sim_homo$Sigma <- diag(c(.01,.02,.05,.04))

s <- sqrt(diag(th_y_sim_homo$Sigma))
th_y_sim_hetero <- th_y_sim_homo
th_y_sim_hetero$paramModel <- 'powLawMixUncorrHetero'
th_y_sim_hetero$kappa <- 0.02

hp_homo <- list(paramModel='powLawMixUncorrHomo')
hp_homo$J <- 2
hp_homo$K <- 2
hp_homo$M <- c(2,2)
hp_hetero <- hp_homo
hp_hetero$paramModel <- 'powLawMixUncorrHetero'

xmin <- 0
xmax <- 20
th_x_sim <- list(fitType='uniform',xmin=xmin,xmax=xmax)

#sim_homo <- simMixedCumProbit(th_y_sim_homo,100,th_x_sim,hp_homo)
sim_homo <- simPowLawMix(th_y_sim_homo,th_x_sim,100,hp_homo)
sim_hetero <- simPowLawMix(th_y_sim_hetero,th_x_sim,100,hp_hetero)


th_y_sim_vect_homo <- theta_y_list2vect(th_y_sim_homo)
th_y_sim_vect_hetero <- theta_y_list2vect(th_y_sim_hetero)

th_y_sim_vect_homo_bar <- theta_y_constr2unconstr(th_y_sim_vect_homo,hp_homo)
th_y_sim_vect_hetero_bar <- theta_y_constr2unconstr(th_y_sim_vect_hetero,hp_hetero)

# Numerically check the gradient calculation
numGrad <- function(th_y,x,Y,hp,transformVar) {
  eps <- 1e-8 # The step size for the finite difference
  f0 <- powLawMixNegLogLik(th_y,x,Y,hp,transformVar)
  N <- length(th_y) # number of variables
  gradVect <- rep(NA,N) # The gradient vector
  # iterate over variables to calculate the numerical gradient
  for(n in 1:N) {
    th_y_eps <- th_y
    th_y_eps[n] <- th_y_eps[n] + eps
    gradVect[n] <- (powLawMixNegLogLik(th_y_eps,x,Y,hp,transformVar) - f0)/eps
  }
  return(gradVect)
}

expect_equal(
  powLawMixGradNegLogLik(th_y_sim_vect_homo,sim_homo$x,sim_homo$Y,hp_homo,F),
  numGrad(th_y_sim_vect_homo,sim_homo$x,sim_homo$Y,hp_homo,F),
  tol=1e-4
)

expect_equal(
  powLawMixGradNegLogLik(th_y_sim_vect_hetero,sim_hetero$x,sim_hetero$Y,hp_hetero,F),
  numGrad(th_y_sim_vect_hetero,sim_hetero$x,sim_hetero$Y,hp_hetero,F),
  tol=1e-4
)

expect_equal(
  powLawMixGradNegLogLik(th_y_sim_vect_homo_bar,sim_homo$x,sim_homo$Y,hp_homo,T),
  numGrad(th_y_sim_vect_homo_bar,sim_homo$x,sim_homo$Y,hp_homo,T),
  tol=1e-4
)

expect_equal(
  powLawMixGradNegLogLik(th_y_sim_vect_hetero_bar,sim_hetero$x,sim_hetero$Y,hp_hetero,T),
  numGrad(th_y_sim_vect_hetero_bar,sim_hetero$x,sim_hetero$Y,hp_hetero,T),
  tol=1e-4
)
