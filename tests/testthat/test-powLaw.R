# Test the functions in powLaw

# th_v has ordering [r,a,b,s,kap]
th_w_homo   <- c(2,.45,1.2,1)
th_w_hetero <- c(2,.45,1.2,1,0.02)
th_w_bar_homo   <- c(log(2),log(.45),1.2,log(1))
th_w_bar_hetero <- c(log(2),log(.45),1.2,log(1),log(0.02))

# Test helper function is_th_w_hetero

expect_error(
  is_th_w_hetero(c(th_w_hetero,10)),
  'Length of th_w should be 4 or 5, not 6'
)

expect_equal(
  is_th_w_hetero(th_w_homo),
  F
)

expect_equal(
  is_th_w_hetero(th_w_hetero),
  T
)

# test powLaw
expect_equal(
  powLaw(5,th_w_homo),
  2*5^.45 + 1.2
)

expect_equal(
  powLaw(5,th_w_hetero),
  2*5^.45 + 1.2
)

# test powLawSigma
expect_equal(
  powLawSigma(5,th_w_homo),
  1
)

expect_equal(
  powLawSigma(5,th_w_hetero,hetSpec='linearSd'),
  1*(1 + 0.02*5)
)

# test powLawDensity
mean_homo <- powLaw(5,th_w_homo)
sig_homo  <- powLawSigma(5,th_w_homo)
mean_hetero <- powLaw(5,th_w_hetero)
sig_hetero  <- powLawSigma(5,th_w_hetero,hetSpec='linearSd')

expect_equal(
  powLawDensity(5,1.22,th_w_homo,hetSpec='none'),
  dnorm(1.22,mean=mean_homo,sd=sig_homo)
)

expect_equal(
  powLawDensity(5,1.22,th_w_hetero,hetSpec='linearSd'),
  dnorm(1.22,mean=mean_hetero,sd=sig_hetero)
)

# test powLawNegLogLikVect
x <- c(1.3,2.1)
w <- c(.7,3.8)
mean_homo <- powLaw(x,th_w_homo)
sig_homo  <- powLawSigma(x,th_w_homo)
mean_hetero <- powLaw(x,th_w_hetero)
sig_hetero  <- powLawSigma(x,th_w_hetero,hetSpec='linearSd')

eta_vect_homo <- powLawNegLogLikVect(th_w_homo,x,w)
expect_equal(
  eta_vect_homo,
  -log(dnorm(w,mean_homo,sig_homo))
)

eta_vect_hetero <- powLawNegLogLikVect(th_w_hetero,x,w,hetSpec='linearSd')
expect_equal(
  eta_vect_hetero,
  -log(dnorm(w,mean_hetero,sig_hetero))
)

# test powLawNegLogLik
expect_equal(
  powLawNegLogLik(th_w_homo,x,w),
  sum(eta_vect_homo)
)

expect_equal(
  powLawNegLogLik(th_w_hetero,x,w,hetSpec='linearSd'),
  sum(eta_vect_hetero)
)

# Numerically check the gradient calculation
numGrad <- function(th_w,x,w,hetSpec,transformVar) {
  eps <- 1e-8 # The step size for the finite difference
  eta0 <- powLawNegLogLik(th_w,x,w,hetSpec,transformVar)
  N <- length(th_w) # number of variables
  gradVect <- rep(NA,N) # The gradient vector
  # iterate over variables to calculate the numerical gradient
  for(n in 1:N) {
    th_w_eps <- th_w
    th_w_eps[n] <- th_w_eps[n] + eps
    gradVect[n] <- (powLawNegLogLik(th_w_eps,x,w,hetSpec,transformVar)-eta0)/eps
  }
  return(gradVect)
}

# homoskedastic / constrained variables
expect_equal(
  powLawGradNegLogLik(th_w_homo,x,w,hetSpec='none',transformVar=F),
  numGrad(th_w_homo,x,w,hetSpec='none',transformVar=F),
  tolerance=1e-4
)

# heteroskedastic / constrained variables
expect_equal(
  powLawGradNegLogLik(th_w_hetero,x,w,hetSpec='linearSd',transformVar=F),
  numGrad(th_w_hetero,x,w,hetSpec='linearSd',transformVar=F),
  tolerance=1e-4
)

# homoskedastic / unconstrained variables
expect_equal(
  powLawGradNegLogLik(th_w_bar_homo,x,w,hetSpec='none',transformVar=T),
  numGrad(th_w_bar_homo,x,w,hetSpec='none',transformVar=T),
  tolerance=1e-4
)

# heteroskedastic / unconstrained variables
expect_equal(
  powLawGradNegLogLik(th_w_bar_hetero,x,w,hetSpec='linearSd',transformVar=T),
  numGrad(th_w_bar_hetero,x,w,hetSpec='linearSd',transformVar=T),
  tolerance=1e-4
)

# test simPowLaw

# A uniform prior on x on the interval 0 to 80
th_x <- c(0,80)
N <- 100

# From random.org between 1 and 1,000,000
set.seed(180190)

# Check simulation for homoskedastic case
expect_error(
  sim_homo <- simPowLaw(N,th_x,th_w_homo,hetSpec='none'),
  NA
)

expect_equal(
  names(sim_homo),
  c('x','w')
)

expect_equal(
  length(sim_homo$x),
  N
)

expect_equal(
  length(sim_homo$w),
  N
)

# Check simulation for heteroskedastic case
expect_error(
  sim_hetero <- simPowLaw(N,th_x,th_w_hetero,hetSpec='linearSd'),
  NA
)

expect_equal(
  names(sim_hetero),
  c('x','w')
)

expect_equal(
  length(sim_hetero$x),
  N
)

expect_equal(
  length(sim_hetero$w),
  N
)

# test fitPowLaw
# Check fit for homoskedastic case
expect_error(
  fit_homo <- fitPowLaw(sim_homo$x,sim_homo$w,hetSpec='none'),
  NA
)

# Check fit for heteroskedastic case
expect_error(
  fit_hetero <- fitPowLaw(sim_hetero$x,sim_hetero$w,hetSpec='linearSd'),
  NA
)