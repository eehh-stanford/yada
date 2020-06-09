# Test the functions in powLaw

# th_v has ordering [r,a,b,s,kap]
th_w_homo   <- c(2,.45,1.2,1.1)
th_w_hetero <- c(2,.45,1.2,1.1,0.02)
th_w_bar_homo   <- c(log(2),log(.45),1.2,log(1.1))
th_w_bar_hetero <- c(log(2),log(.45),1.2,log(1.1),log(0.02))

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
x <- c(1.3,2.1)
w <- c(.7,3.8)

expect_equal(
  powLaw(x,th_w_homo,transformVar=F),
  2*x^.45 + 1.2
)

expect_equal(
  powLaw(x,th_w_bar_homo,transformVar=T),
  2*x^.45 + 1.2
)

expect_equal(
  powLaw(x,th_w_hetero,transformVar=F),
  2*x^.45 + 1.2
)

expect_equal(
  powLaw(x,th_w_bar_hetero,transformVar=T),
  2*x^.45 + 1.2
)

# test powLawSigma
expect_equal(
  powLawSigma(x,th_w_homo,hetSpec='none',transformVar=F),
  1.1 # a scalar is expected, even though x is a vector
)

expect_equal(
  powLawSigma(x,th_w_bar_homo,hetSpec='none',transformVar=T),
  1.1 # a scalar is expected, even though x is a vector
)

expect_equal(
  powLawSigma(x,th_w_hetero,hetSpec='sd_x',transformVar=F),
  1.1*(1 + 0.02*x)
)

expect_equal(
  powLawSigma(x,th_w_bar_hetero,hetSpec='sd_x',transformVar=T),
  1.1*(1 + 0.02*x)
)

expect_equal(
  powLawSigma(x,th_w_hetero,hetSpec='sd_resp',transformVar=F),
  1.1*(1 + 0.02*2*x^.45)
)

expect_equal(
  powLawSigma(x,th_w_bar_hetero,hetSpec='sd_resp',transformVar=T),
  1.1*(1 + 0.02*2*x^.45)
)

# test powLawDensity
mean_homo <- powLaw(5,th_w_homo)
sig_homo  <- powLawSigma(5,th_w_homo)
mean_hetero <- powLaw(5,th_w_hetero)
sig_heterox  <- powLawSigma(5,th_w_hetero,hetSpec='sd_x')
sig_heteror  <- powLawSigma(5,th_w_hetero,hetSpec='sd_resp')

expect_equal(
  powLawDensity(5,1.22,th_w_homo,hetSpec='none'),
  dnorm(1.22,mean=mean_homo,sd=sig_homo)
)

expect_equal(
  powLawDensity(5,1.22,th_w_hetero,hetSpec='sd_x'),
  dnorm(1.22,mean=mean_hetero,sd=sig_heterox)
)

expect_equal(
  powLawDensity(5,1.22,th_w_hetero,hetSpec='sd_resp'),
  dnorm(1.22,mean=mean_hetero,sd=sig_heteror)
)

# Directly calculate the log likelihood, homoskedastic case
mean_homo <- powLaw(x,th_w_homo)
sig_homo  <- powLawSigma(x,th_w_homo)

# -- with transformVar = F
eta_vect_homo <- powLawNegLogLikVect(th_w_homo,x,w,hetSpec='none',transformVar=F)
expect_equal(
  eta_vect_homo,
  -log(dnorm(w,mean_homo,sig_homo))
)

expect_equal(
  powLawNegLogLik(th_w_homo,x,w,hetSpec='none',transformVar=F),
  sum(eta_vect_homo)
)

# -- with transformVar = T
eta_vect_homo <- powLawNegLogLikVect(th_w_bar_homo,x,w,hetSpec='none',transformVar=T)
expect_equal(
  eta_vect_homo,
  -log(dnorm(w,mean_homo,sig_homo))
)

# Directly calculate the log likelihood, hetSpec = 'sd_x'
mean_hetero <- powLaw(x,th_w_hetero)
sig_heterox  <- powLawSigma(x,th_w_hetero,hetSpec='sd_x')

# -- with transformVar = F
eta_vect_heterox <- powLawNegLogLikVect(th_w_hetero,x,w,hetSpec='sd_x',transformVar=F)
expect_equal(
  eta_vect_heterox,
  -log(dnorm(w,mean_hetero,sig_heterox))
)

expect_equal(
  powLawNegLogLik(th_w_hetero,x,w,hetSpec='sd_x'),
  sum(eta_vect_heterox)
)

# -- with transformVar = T
eta_vect_heterox <- powLawNegLogLikVect(th_w_bar_hetero,x,w,hetSpec='sd_x',transformVar=T)
expect_equal(
  eta_vect_heterox,
  -log(dnorm(w,mean_hetero,sig_heterox))
)

expect_equal(
  powLawNegLogLik(th_w_hetero,x,w,hetSpec='sd_x'),
  sum(eta_vect_heterox)
)

# Directly calculate the log likelihood, hetSpec = 'sd_resp'
sig_heteror  <- powLawSigma(x,th_w_hetero,hetSpec='sd_resp')

# -- with transformVar = F
eta_vect_heteror <- powLawNegLogLikVect(th_w_hetero,x,w,hetSpec='sd_resp',transformVar=F)
expect_equal(
  eta_vect_heteror,
  -log(dnorm(w,mean_hetero,sig_heteror))
)

expect_equal(
  powLawNegLogLik(th_w_hetero,x,w,hetSpec='sd_resp',transformVar=F),
  sum(eta_vect_heteror)
)

# -- with transformVar = T
eta_vect_heteror <- powLawNegLogLikVect(th_w_bar_hetero,x,w,hetSpec='sd_resp',transformVar=T)
expect_equal(
  eta_vect_heteror,
  -log(dnorm(w,mean_hetero,sig_heteror))
)

expect_equal(
  powLawNegLogLik(th_w_bar_hetero,x,w,hetSpec='sd_resp',transformVar=T),
  sum(eta_vect_heteror)
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

# Check simulation for heteroskedastic case (sd_x)
expect_error(
  sim_heterox <- simPowLaw(N,th_x,th_w_hetero,hetSpec='sd_x'),
  NA
)

expect_equal(
  names(sim_heterox),
  c('x','w')
)

expect_equal(
  length(sim_heterox$x),
  N
)

expect_equal(
  length(sim_heterox$w),
  N
)

# Check simulation for heteroskedastic case (sd_resp)
expect_error(
  sim_heteror <- simPowLaw(N,th_x,th_w_hetero,hetSpec='sd_resp'),
  NA
)

expect_equal(
  names(sim_heteror),
  c('x','w')
)

expect_equal(
  length(sim_heteror$x),
  N
)

expect_equal(
  length(sim_heteror$w),
  N
)

# test fitPowLaw
# Check fit for homoskedastic case
expect_error(
  fit_homo <- fitPowLaw(sim_homo$x,sim_homo$w,hetSpec='none'),
  NA
)

# Check fit for heteroskedastic case (sd_x)
expect_error(
  fit_heterox <- fitPowLaw(sim_heterox$x,sim_heterox$w,hetSpec='sd_x'),
  NA
)

# Check fit for heteroskedastic case (sd_resp)
expect_error(
  fit_heteror <- fitPowLaw(sim_heteror$x,sim_heteror$w,hetSpec='sd_resp'),
  NA
)
