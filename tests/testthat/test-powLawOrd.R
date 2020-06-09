# Test the functions in powLawOrd

# th_v has ordering [rho,tau_1,...tau_2,s,kap]
th_v_homo   <- c(.75,15,0.01)
th_v_hetero <- c(.75,15,0.01,0.02)
th_v_bar_homo   <- c(log(.75),15,log(0.01))
th_v_bar_hetero <- c(log(.75),15,log(0.01),log(0.02))

# test powLawOrd
x0 <- 10
x1 <- 20
x2 <- 30
x <- c(x0,x1,x2)
v <- c(0,1,2)

expect_equal(
  powLawOrd(x,th_v_homo,transformVar=F),
  x^.75
)

expect_equal(
  powLawOrd(x,th_v_bar_homo,transformVar=T),
  x^.75
)

expect_equal(
  powLawOrd(x,th_v_hetero,transformVar=F),
  x^.75
)

expect_equal(
  powLawOrd(x,th_v_bar_hetero,transformVar=T),
  x^.75
)

# test powLawOrdSigma
expect_equal(
  powLawOrdSigma(x,th_v_homo,hetSpec='none',transformVar=F),
  0.01
)

expect_equal(
  powLawOrdSigma(x,th_v_bar_homo,hetSpec='none',transformVar=T),
  0.01
)

expect_equal(
  powLawOrdSigma(x,th_v_hetero,hetSpec='sd_x',transformVar=F),
  0.01*(1+x*0.02)
)

expect_equal(
  powLawOrdSigma(x,th_v_bar_hetero,hetSpec='sd_x',transformVar=T),
  0.01*(1+x*0.02)
)

expect_equal(
  powLawOrdSigma(x,th_v_hetero,hetSpec='sd_resp',transformVar=F),
  0.01*(1+x^.75*0.02)
)

expect_equal(
  powLawOrdSigma(x,th_v_bar_hetero,hetSpec='sd_resp',transformVar=T),
  0.01*(1+x^.75*0.02)
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


# Directly calculate the log likelihood, homoskedastic case
eta_v_homo0 <- -log(pnorm( (tau1 - x0^rho)/s ))
eta_v_homo1 <- -log( pnorm( (tau2 - x1^rho)/s ) - pnorm( (tau1 - x1^rho)/s ))
eta_v_homo2 <- -log( 1 - pnorm( (tau2 - x2^rho)/s ))

# -- with transformVar = F
eta_v_vect_homo <- powLawOrdNegLogLikVect(th_v_homo2,x,v,hetSpec='none',transformVar=F)
expect_equal(
  eta_v_vect_homo,
 c(eta_v_homo0,eta_v_homo1,eta_v_homo2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_homo2,x,v,hetSpec='none',transformVar=F),
  sum(eta_v_vect_homo)
)

# -- with transformVar = T
eta_v_vect_homo <- powLawOrdNegLogLikVect(th_v_bar_homo2,x,v,hetSpec='none',transformVar=T)
expect_equal(
  eta_v_vect_homo,
 c(eta_v_homo0,eta_v_homo1,eta_v_homo2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_bar_homo2,x,v,hetSpec='none',transformVar=T),
  sum(eta_v_vect_homo)
)

# Directly calculate the log likelihood, hetSpec = 'sd_x'
eta_v_hetero0 <- -log(pnorm( (tau1 - x0^rho)/(s*(1+kap*x0)) ))
eta_v_hetero1 <- -log( pnorm( (tau2 - x1^rho)/(s*(1+kap*x1)) ) - pnorm( (tau1 - x1^rho)/(s*(1+kap*x1)) ))
eta_v_hetero2 <- -log( 1 - pnorm( (tau2 - x2^rho)/(s*(1+kap*x2)) ))

# -- with transformVar = F
eta_v_vect_hetero <- powLawOrdNegLogLikVect(th_v_hetero2,x,v,hetSpec='sd_x',transformVar=F)
expect_equal(
  eta_v_vect_hetero,
 c(eta_v_hetero0,eta_v_hetero1,eta_v_hetero2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_hetero2,x,v,hetSpec='sd_x',transformVar=F),
  sum(eta_v_vect_hetero)
)

# -- with transformVar = T
eta_v_vect_hetero <- powLawOrdNegLogLikVect(th_v_bar_hetero2,x,v,hetSpec='sd_x',transformVar=T)
expect_equal(
  eta_v_vect_hetero,
 c(eta_v_hetero0,eta_v_hetero1,eta_v_hetero2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_bar_hetero2,x,v,hetSpec='sd_x',transformVar=T),
  sum(eta_v_vect_hetero)
)

# Directly calculate the log likelihood, hetSpec = 'sd_resp'
eta_v_hetero0 <- -log(pnorm( (tau1 - x0^rho)/(s*(1+kap*x0^rho)) ))
eta_v_hetero1 <- -log( pnorm( (tau2 - x1^rho)/(s*(1+kap*x1^rho)) ) - pnorm( (tau1 - x1^rho)/(s*(1+kap*x1^rho)) ))
eta_v_hetero2 <- -log( 1 - pnorm( (tau2 - x2^rho)/(s*(1+kap*x2^rho)) ))

# -- with transformVar = F
eta_v_vect_hetero <- powLawOrdNegLogLikVect(th_v_hetero2,x,v,hetSpec='sd_resp',transformVar=F)
expect_equal(
  eta_v_vect_hetero,
 c(eta_v_hetero0,eta_v_hetero1,eta_v_hetero2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_hetero2,x,v,hetSpec='sd_resp',transformVar=F),
  sum(eta_v_vect_hetero)
)

# -- with transformVar = T
eta_v_vect_hetero <- powLawOrdNegLogLikVect(th_v_bar_hetero2,x,v,hetSpec='sd_resp',transformVar=T)
expect_equal(
  eta_v_vect_hetero,
 c(eta_v_hetero0,eta_v_hetero1,eta_v_hetero2) 
)

expect_equal(
  powLawOrdNegLogLik(th_v_bar_hetero2,x,v,hetSpec='sd_resp',transformVar=T),
  sum(eta_v_vect_hetero)
)

# test calc_M
expect_equal(
  calc_M(th_v_homo,hetSpec='none'),
  1
)

expect_equal(
  calc_M(th_v_hetero,hetSpec='sd_x'),
  1
)

expect_equal(
  calc_M(th_v_hetero,hetSpec='sd_resp'),
  1
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

# Check simulation for heteroskedastic case (sd_x)
expect_error(
  sim_heterox <- simPowLawOrd(N,th_x,th_v_hetero,hetSpec='sd_x'),
  NA
)

expect_equal(
  names(sim_heterox),
  c('x','v','vstar')
)

expect_equal(
  length(sim_heterox$x),
  N
)

expect_equal(
  length(sim_heterox$v),
  N
)

expect_equal(
  length(sim_heterox$vstar),
  N
)

# Check simulation for heteroskedastic case (sd_resp)
expect_error(
  sim_heteror <- simPowLawOrd(N,th_x,th_v_hetero,hetSpec='sd_resp'),
  NA
)

expect_equal(
  names(sim_heteror),
  c('x','v','vstar')
)

expect_equal(
  length(sim_heteror$x),
  N
)

expect_equal(
  length(sim_heteror$v),
  N
)

expect_equal(
  length(sim_heteror$vstar),
  N
)

# test fitPowLawOrd
# Check fit for homoskedastic case
expect_error(
  fit_homo <- fitPowLawOrd(sim_homo$x,sim_homo$v,hetSpec='none'),
  NA
)

# Check fit for heteroskedastic case (sd_x)
expect_error(
  fit_heterox <- fitPowLawOrd(sim_heterox$x,sim_heterox$v,hetSpec='sd_x'),
  NA
)

# Check fit for heteroskedastic case (sd_resp)
expect_error(
  fit_heteror <- fitPowLawOrd(sim_heteror$x,sim_heteror$v,hetSpec='sd_resp'),
  NA
)
