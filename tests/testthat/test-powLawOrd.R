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
  powLawOrdSigma(5,th_v_homo,F,F),
  0.01
)

expect_equal(
  powLawOrdSigma(5,th_v_bar_homo,F,T),
  0.01
)

expect_equal(
  powLawOrdSigma(5,th_v_hetero,T,F),
  0.01*(1+5*0.02)
)

expect_equal(
  powLawOrdSigma(5,th_v_bar_hetero,T,T),
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
hp_homo <- list(paramModel='powLawOrdHomo')
hp_homo$J <- 1
hp_homo$M <- length(tau)

hp_hetero <- list(paramModel='powLawOrdHetero')
hp_hetero$J <- 1
hp_hetero$M <- length(tau)

x0 <- 10
x1 <- 20
x2 <- 30
x_list <- list()
x_list[[1]] <- x0
x_list[[2]] <- x1
x_list[[3]] <- x2

# Directly calculate the log likelihood homoskedastic case
eta_v_homo <- -log(pnorm( (tau1 - x0^rho)/s ))
eta_v_homo <- eta_v_homo - log( pnorm( (tau2 - x1^rho)/s ) - pnorm( (tau1 - x1^rho)/s ))
eta_v_homo <- eta_v_homo - log( 1 - pnorm( (tau2 - x2^rho)/s ))

expect_equal(
  powLawOrdNegLogLik(th_v_homo2,x_list,hetero=F,transformVar=F),
  eta_v_homo
)

expect_equal(
  powLawOrdNegLogLik(th_v_bar_homo2,x_list,hetero=F,transformVar=T,hp=hp_homo),
  eta_v_homo
)

# Directly calculate the log likelihood heteroskedastic case
eta_v_hetero <- -log(pnorm( (tau1 - x0^rho)/(s*(1+kap*x0)) ))
eta_v_hetero <- eta_v_hetero - log( pnorm( (tau2 - x1^rho)/(s*(1+kap*x1)) ) - pnorm( (tau1 - x1^rho)/(s*(1+kap*x1)) ))
eta_v_hetero <- eta_v_hetero - log( 1 - pnorm( (tau2 - x2^rho)/(s*(1+kap*x2)) ))

expect_equal(
  powLawOrdNegLogLik(th_v_hetero2,x_list,hetero=T,transformVar=F),
  eta_v_hetero
)

expect_equal(
  powLawOrdNegLogLik(th_v_bar_hetero2,x_list,hetero=T,transformVar=T,hp=hp_hetero),
  eta_v_hetero
)


# test simPowLawOrd
N <- 100 # Number of points for simulation

# A uniform prior on x on the interval 0 to 80
th_x <- c(0,80)

# Check simulation for homoskedastic case
expect_error(
  sim_homo <- simPowLawOrd(N,th_x,th_v_homo,hetero=F),
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
  sim_hetero <- simPowLawOrd(N,th_x,th_v_hetero,hetero=T),
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
  fit_homo <- fitPowLawOrd(sim_homo$x,sim_homo$v,hetero=F),
  NA
)

# Check fit for heteroskedastic case
expect_error(
  fit_hetero <- fitPowLawOrd(sim_homo$x,sim_homo$v,hetero=T),
  NA
)
