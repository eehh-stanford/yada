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

# Parameters for second ordinal variable
rho2 <- .7
tau2 <- c(1.3,3.1)
M2 <- 2
s2 <- .4
kappa2 <- .015

# Parameters for first continuous variable
a1 <- 10
r1 <- .6
b1 <- 5
s3 <- 3
kappa3 <- .01

# Parameters for second continuous variable
a2 <- 8
r2 <- .4
b2 <- -2
s4 <- 1.5
kappa4 <- .02

N <- 100
th_x <- c(0,20)

# Create simulated data for all four variables and various choices for the
# heteroskedasticity.
sim1_none <- simPowLawOrd(N,th_x,c(rho1,tau1,s1),       hetSpec='none')
sim1_sdx  <- simPowLawOrd(N,th_x,c(rho1,tau1,s1,kappa1),hetSpec='sd_x')
sim1_sdr  <- simPowLawOrd(N,th_x,c(rho1,tau1,s1,kappa1),hetSpec='sd_resp')

sim2_none <- simPowLawOrd(N,th_x,c(rho2,tau2,s2),       hetSpec='none')
sim2_sdx  <- simPowLawOrd(N,th_x,c(rho2,tau2,s2,kappa2),hetSpec='sd_x')
sim2_sdr  <- simPowLawOrd(N,th_x,c(rho2,tau2,s2,kappa2),hetSpec='sd_resp')

sim3_none <- simPowLaw(N,th_x,c(a1,r1,b1,s3),       hetSpec='none')
sim3_sdx  <- simPowLaw(N,th_x,c(a1,r1,b1,s3,kappa3),hetSpec='sd_x')
sim3_sdr  <- simPowLaw(N,th_x,c(a1,r1,b1,s3,kappa3),hetSpec='sd_resp')

sim4_none <- simPowLaw(N,th_x,c(a2,r2,b2,s4),       hetSpec='none')
sim4_sdx  <- simPowLaw(N,th_x,c(a2,r2,b2,s4,kappa4),hetSpec='sd_x')
sim4_sdr  <- simPowLaw(N,th_x,c(a2,r2,b2,s4,kappa4),hetSpec='sd_resp')

# Single variable, ordinal, hetSpec = 'none'
modSpec1_none <- list(meanSpec='powLaw')
modSpec1_none$J <- 1
modSpec1_none$M <- M1
modSpec1_none$hetSpec <- 'none'
expect_error(
  fit <- fit_theta_y(sim1_none$x,sim1_none$v,modSpec1_none),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_x'
modSpec1_sdx <- list(meanSpec='powLaw')
modSpec1_sdx$J <- 1
modSpec1_sdx$M <- M1
modSpec1_sdx$hetSpec <- 'sd_x'
modSpec1_sdx$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(sim1_sdx$x,sim1_sdx$v,modSpec1_sdx),
  NA
)

# Single variable, ordinal, hetSpec = 'sd_resp'
modSpec1_sdr <- list(meanSpec='powLaw')
modSpec1_sdr$J <- 1
modSpec1_sdr$M <- M1
modSpec1_sdr$hetSpec <- 'sd_resp'
modSpec1_sdr$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(sim1_sdr$x,sim1_sdr$v,modSpec1_sdr),
  NA
)

# Single variable, continuous, hetSpec = 'none'
modSpec3_none <- list(meanSpec='powLaw')
modSpec3_none$K <- 1
modSpec3_none$hetSpec <- 'none'
expect_error(
  fit <- fit_theta_y(sim3_none$x,sim3_none$w,modSpec3_none),
  NA
)

# Single variable, continuous, hetSpec = 'sd_x'
modSpec3_sdx <- list(meanSpec='powLaw')
modSpec3_sdx$K <- 1
modSpec3_sdx$hetSpec <- 'sd_x'
modSpec3_sdx$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(sim3_sdx$x,sim3_sdx$w,modSpec3_sdx),
  NA
)

# Single variable, continuous, hetSpec = 'sd_resp'
modSpec3_sdr <- list(meanSpec='powLaw')
modSpec3_sdr$K <- 1
modSpec3_sdr$hetSpec <- 'sd_resp'
modSpec3_sdr$hetGroups <- 1
expect_error(
  fit <- fit_theta_y(sim3_sdr$x,sim3_sdr$w,modSpec3_sdr),
  NA
)


## Four variables, homoskedastic
#modSpec <- list(meanSpec='powLaw')
#modSpec$J <- 2
#modSpec$K <- 2
#modSpec$hetSpec <- 'none'
#modSpec$cdepSpec <- 'indep'
#
#expect_error(
#  fit <- fit_theta_y(sim3_sdr$x,sim3_sdr$w,modSpec3_sdr),
#  NA
#)

