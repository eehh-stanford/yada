# Test fit_theta_y

# From random.org between 1 and 1,000,000:
set.seed(874287)

N <- 100

# Parameters for first ordinal variable
rho1 <- .5
tau1 <- c(1.5,3.5)
M1 <- 2
s1 <- .5
kappa1 <- .01

# Parameters for first continuous variable
a1 <- 10
r1 <- .6
b1 <- 5
s3 <- 3
kappa3 <- .01


N <- 100
th_x <- c(0,20)

# Single variable, ordinal, hetSpec = 'none'
th_v <- c(rho1,tau1,s1)
sim <- simPowLawOrd(N,th_x,th_v)
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- M1
modSpec$hetSpec <- 'none'

expect_error(
  fit <- fit_theta_y(sim$x,sim$v,modSpec),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_x'
th_v <- c(rho1,tau1,s1,kappa1)
sim <- simPowLawOrd(N,th_x,th_v,hetSpec='sd_x')
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- M1
modSpec$hetSpec <- 'sd_x'
modSpec$hetGroups <- 1

expect_error(
  fit <- fit_theta_y(sim$x,sim$v,modSpec),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_resp'
sim <- simPowLawOrd(N,th_x,th_v,hetSpec='sd_x')
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- M1
modSpec$hetSpec <- 'sd_resp'
modSpec$hetGroups <- 1

expect_error(
  fit <- fit_theta_y(sim$x,sim$v,modSpec),
  NA
)

# Single variable, continuous, hetSpec = 'none'
th_w <- c(a1,r1,b1,s3)
sim <- simPowLaw(N,th_x,th_w)
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'none'

expect_error(
  fit <- fit_theta_y(sim$x,sim$w,modSpec),
  NA
)

# Single variable, continuous, hetSpec = 'sd_x'
th_w <- c(a1,r1,b1,s3,kappa3)
sim <- simPowLaw(N,th_x,th_w,hetSpec='sd_x')
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'sd_x'
modSpec$hetGroups <- 1

expect_error(
  fit <- fit_theta_y(sim$x,sim$w,modSpec),
  NA
)

# Single variable, continuous, hetSpec = 'sd_x'
th_w <- c(a1,r1,b1,s3,kappa3)
sim <- simPowLaw(N,th_x,th_w,hetSpec='sd_resp')
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'sd_resp'
modSpec$hetGroups <- 1

expect_error(
  fit <- fit_theta_y(sim$x,sim$w,modSpec),
  NA
)


if(F) {
x <- seq(0,10,len=6)
v1 <- c(0,0,0,1,1,1)
v2 <- c(0,0,1,1,1,1)
# w1 made with: 10*x^.5 + 5 + rnorm(length(x))
w1 <- c(4.701958,19.452807,24.210377,28.544634,32.488103,36.337787)
# w2 made with: 5*x^.75 - 1 + rnorm(length(x))/2
w2 <- c(0.2570962,7.8094604,13.1529883,17.6184679,22.9311745,26.9613796)


# Single variable, ordinal, hetSpec = 'none'
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- 1
modSpec$hetSpec <- 'none'
expect_error(
  fit <- fit_theta_y(x,v1,modSpec),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_x'
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- 1
modSpec$hetSpec <- 'sd_x'
modSpec$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(x,v1,modSpec),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_resp'
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- 1
modSpec$hetSpec <- 'sd_resp'
modSpec$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(x,v1,modSpec),
  NA
)


# Single variable, continuous, hetSpec = 'sd_x'
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'sd_x'
modSpec$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(x,w1,modSpec),
  NA
)

# Single variable, continuous, hetSpec = 'sd_resp'
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'sd_resp'
modSpec$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(x,w1,modSpec),
  NA
)
}
