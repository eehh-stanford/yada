# Test the functions in powLawOrd

# th_v has ordering [rho,tau_1,...tau_2,s,kap]
th_v_homo   <- c(.75,15,0.01)
th_v_hetero <- c(.75,15,0.01,0.02)
th_v_bar_homo   <- c(log(.75),15,log(0.01))
th_v_bar_hetero <- c(log(.75),15,log(0.01),log(0.02))

# test powLawOrd
expect_equal(
  powLawOrd(5,th_v_homo),
  5^.75
)

expect_equal(
  powLawOrd(5,th_v_bar_homo,T),
  5^.75
)

expect_equal(
  powLawOrd(5,th_v_hetero),
  5^.75
)

expect_equal(
  powLawOrd(5,th_v_bar_hetero,T),
  5^.75
)

# test powLawOrdSigma
expect_equal(
  powLawOrdSigma(5,th_v_homo,'none',F),
  0.01
)

expect_equal(
  powLawOrdSigma(5,th_v_bar_homo,'none',T),
  0.01
)

expect_equal(
  powLawOrdSigma(5,th_v_hetero,'linearSd',F),
  0.01*(1+5*0.02)
)

expect_equal(
  powLawOrdSigma(5,th_v_bar_hetero,'linearSd',T),
  0.01*(1+5*0.02)
)

# test powLawOrdNegLogLik

rho <- 0.90
tau1 <- 10
tau2 <- 14
tau <- c(tau1,tau2)
  s <- 1.9
kap <- 0.01

th_v_homo2 <- c(rho,tau,s)
th_v_hetero2 <- c(rho,tau,s,kap)
th_v_bar_homo2 <- c(log(rho),tau[1],log(tau[2]-tau[1]),log(s))
th_v_bar_hetero2 <- c(log(rho),tau[1],log(tau[2]-tau[1]),log(s),log(kap))

x0 <- 10
x1 <- 20
x2 <- 30
x <- c(x0,x1,x2)
v <- c(0,1,2)

# Directly calculate the log likelihood, homoskedastic case
eta_v_homo0 <- -log(pnorm( (tau1 - x0^rho)/s ))
eta_v_homo1 <- -log( pnorm( (tau2 - x1^rho)/s ) - pnorm( (tau1 - x1^rho)/s ))
eta_v_homo2 <- -log( 1 - pnorm( (tau2 - x2^rho)/s ))

eta_v_vect_homo <- powLawOrdNegLogLikVect(th_v_homo2,x,v,hetSpec='none',transformVar=F)
expect_equal(
  eta_v_vect_homo,
 c(eta_v_homo0,eta_v_homo1,eta_v_homo2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_homo2,x,v,hetSpec='none',transformVar=F),
  sum(eta_v_vect_homo)
)

# Directly calculate the log likelihood, heteroskedastic case
eta_v_hetero0 <- -log(pnorm( (tau1 - x0^rho)/(s*(1+kap*x0)) ))
eta_v_hetero1 <- -log( pnorm( (tau2 - x1^rho)/(s*(1+kap*x1)) ) - pnorm( (tau1 - x1^rho)/(s*(1+kap*x1)) ))
eta_v_hetero2 <- -log( 1 - pnorm( (tau2 - x2^rho)/(s*(1+kap*x2)) ))

eta_v_vect_hetero <- powLawOrdNegLogLikVect(th_v_hetero2,x,v,hetSpec='linearSd',transformVar=F)
expect_equal(
  eta_v_vect_hetero,
 c(eta_v_hetero0,eta_v_hetero1,eta_v_hetero2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_hetero2,x,v,hetSpec='linearSd',transformVar=F),
  sum(eta_v_vect_hetero)
)

# Numerically check the gradient calculation
numGrad <- function(th_v,x,v,hetSpec,transformVar) {
  eps <- 1e-8 # The step size for the finite difference
  eta0 <- powLawOrdNegLogLik(th_v,x,v,hetSpec,transformVar)
  N <- length(th_v) # number of variables
  gradVect <- rep(NA,N) # The gradient vector
  # iterate over variables to calculate the numerical gradient
  for(n in 1:N) {
    th_v_eps <- th_v
    th_v_eps[n] <- th_v_eps[n] + eps
    gradVect[n] <- (powLawOrdNegLogLik(th_v_eps,x,v,hetSpec,transformVar)-eta0)/eps
  }
  return(gradVect)
}

# homoskedastic / constrained variables

expect_equal(
  powLawOrdGradNegLogLik(th_v_homo2,x,v,hetSpec='none',transformVar=F),
  numGrad(th_v_homo2,x,v,hetSpec='none',transformVar=F),
  tolerance=1e-4
)

# heteroskedastic / constrained variables
expect_equal(
  powLawOrdGradNegLogLik(th_v_hetero2,x,v,hetSpec='linearSd',transformVar=F),
  numGrad(th_v_hetero2,x,v,hetSpec='linearSd',transformVar=F),
  tolerance=1e-4
)

# homoskedastic / unconstrained variables
expect_equal(
  powLawOrdGradNegLogLik(th_v_bar_homo2,x,v,hetSpec='none',transformVar=T),
  numGrad(th_v_bar_homo2,x,v,hetSpec='none',transformVar=T),
  tolerance=1e-4
)

# heteroskedastic / unconstrained variables
expect_equal(
  powLawOrdGradNegLogLik(th_v_bar_hetero2,x,v,hetSpec='linearSd',transformVar=T),
  numGrad(th_v_bar_hetero2,x,v,hetSpec='linearSd',transformVar=T),
  tolerance=1e-4
)

# Check the gradient directly for when hetero=F and transformVar=T
delta_Phi0 <- pnorm((tau1-x0^rho)/s)
delta_Phi1 <- pnorm((tau2-x1^rho)/s) - pnorm((tau1-x1^rho)/s)
delta_Phi2 <- 1 - pnorm((tau2-x2^rho)/s)
delta_phi0 <- dnorm((tau1-x0^rho)/s)
delta_phi1 <- dnorm((tau2-x1^rho)/s) - dnorm((tau1-x1^rho)/s)
delta_phi2 <- -dnorm((tau2-x2^rho)/s)

# These are actually the barred variables
grad_rho <- delta_phi0 / delta_Phi0 * x0^rho * log(x0) / s * rho
grad_rho <- grad_rho + delta_phi1 / delta_Phi1 * x1^rho * log(x1) / s * rho
grad_rho <- grad_rho + delta_phi2 / delta_Phi2 * x2^rho * log(x2) / s * rho

grad_tau1 <- - delta_phi0 / delta_Phi0 / s
grad_tau1 <- grad_tau1 - delta_phi1 / delta_Phi1 / s
grad_tau1 <- grad_tau1 - delta_phi2 / delta_Phi2 / s

grad_tau2 <- -dnorm((tau2-x1^rho)/s) / delta_Phi1 / s * (tau2-tau1)
grad_tau2 <- grad_tau2 + dnorm( (tau2-x2^rho)/s) / delta_Phi2 / s * (tau2-tau1)

grad_s <- dnorm((tau1-x0^rho)/s)*(tau1-x0^rho)/delta_Phi0/s
grad_s <- grad_s + ( dnorm((tau2-x1^rho)/s)*(tau2-x1^rho) - dnorm((tau1-x1^rho)/s)*(tau1-x1^rho) ) / delta_Phi1 / s
grad_s <- grad_s - dnorm((tau2-x2^rho)/s)*(tau2-x2^rho) / delta_Phi2 / s

expect_equal(
  powLawOrdGradNegLogLik(th_v_bar_homo2,x,v,hetSpec='none',transformVar=T),
  c(grad_rho,grad_tau1,grad_tau2,grad_s)
)

# test simPowLawOrd
N <- 100 # Number of points for simulation

# A uniform prior on x on the interval 0 to 80
th_x <- c(0,80)

# Check simulation for homoskedastic case
expect_error(
  sim_homo <- simPowLawOrd(N,th_x,th_v_homo,hetSpec='none'),
  NA
)

expect_equal(
  names(sim_homo),
  c('x','v','vstar')
)

expect_equal(
  length(sim_homo$x),
  N
)

expect_equal(
  length(sim_homo$v),
  N
)

expect_equal(
  length(sim_homo$vstar),
  N
)

# Check simulation for heteroskedastic case
expect_error(
  sim_hetero <- simPowLawOrd(N,th_x,th_v_hetero,hetSpec='linearSd'),
  NA
)

expect_equal(
  names(sim_hetero),
  c('x','v','vstar')
)

expect_equal(
  length(sim_hetero$x),
  N
)

expect_equal(
  length(sim_hetero$v),
  N
)

expect_equal(
  length(sim_hetero$vstar),
  N
)

# test fitPowLawOrd
# Check fit for homoskedastic case
expect_error(
  fit_homo <- fitPowLawOrd(sim_homo$x,sim_homo$v,hetSpec='none'),
  NA
)

# Check fit for heteroskedastic case
expect_error(
  fit_hetero <- fitPowLawOrd(sim_hetero$x,sim_hetero$v,hetSpec='linearSd'),
  NA
)
