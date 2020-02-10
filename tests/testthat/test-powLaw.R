# Test the functions in powLawOrd

# th_v has ordering [rho,tau_1,...tau_2,s,kap]
th_w_homo   <- c(2,.45,1.2,1)
th_w_hetero <- c(2,.45,1.2,1,0.02)
th_w_bar_homo   <- c(log(2),log(.45),1.2,log(1))
th_w_bar_hetero <- c(log(2),log(.45),1.2,log(1),log(0.02))

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
  powLawSigma(5,th_w_homo,hetero=F),
  1
)

expect_equal(
  powLawSigma(5,th_w_hetero,hetero=T),
  1*(1 + 0.02*5)
)

# test powLawDensity
mean_homo <- powLaw(5,th_w_homo)
sig_homo  <- powLawSigma(5,th_w_homo,hetero=F)
mean_hetero <- powLaw(5,th_w_hetero)
sig_hetero  <- powLawSigma(5,th_w_hetero,hetero=T)

expect_equal(
  powLawDensity(5,1.22,th_w_homo,hetero=F),
  dnorm(1.22,mean=mean_homo,sd=sig_homo)
)

expect_equal(
  powLawDensity(5,1.22,th_w_hetero,hetero=T),
  dnorm(1.22,mean=mean_hetero,sd=sig_hetero)
)

# test powLawNegLogLik
x <- c(1.3,2.1)
w <- c(.7,3.8)

mean_homo <- powLaw(x,th_w_homo)
sig_homo  <- powLawSigma(x,th_w_homo,hetero=F)
mean_hetero <- powLaw(x,th_w_hetero)
sig_hetero  <- powLawSigma(x,th_w_hetero,hetero=T)

expect_equal(
  powLawNegLogLik(th_w_homo,x,w,hetero=F),
  -sum(log(dnorm(w,mean_homo,sig_homo)))
)

expect_equal(
  powLawNegLogLik(th_w_hetero,x,w,hetero=T),
  -sum(log(dnorm(w,mean_hetero,sig_hetero)))
)

# test simPowLaw

# A uniform prior on x on the interval 0 to 80
th_x <- c(0,80)
N <- 100

# From random.org between 1 and 1,000,000
set.seed(180190)

# Check simulation for homoskedastic case
expect_error(
  sim_homo <- simPowLaw(N,th_x,th_w_homo,hetero=F),
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
  sim_hetero <- simPowLaw(N,th_x,th_w_hetero,hetero=T),
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
  fit_homo <- fitPowLaw(sim_homo$x,sim_homo$w,hetero=F),
  NA
)

# Check fit for heteroskedastic case
expect_error(
  fit_hetero <- fitPowLaw(sim_hetero$x,sim_hetero$w,hetero=T),
  NA
)
