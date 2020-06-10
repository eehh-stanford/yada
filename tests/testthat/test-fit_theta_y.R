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

# Multi-variable, homoskedastic
modSpec_homo <- list(meanSpec='powLaw')
modSpec_homo$J <- 2
modSpec_homo$K <- 2
modSpec_homo$M <- c(2,2)
modSpec_homo$hetSpec <- 'none'
modSpec_homo$cdepSpec <- 'indep'

th_y_sim_homo <- list(modSpec=modSpec_homo)
th_y_sim_homo$rho <- c(.75,.5)
th_y_sim_homo$tau <- list()
th_y_sim_homo$tau[[1]] <- c(2,5)
th_y_sim_homo$tau[[2]] <- c(1,3)
th_y_sim_homo$a <- c(2,3)
th_y_sim_homo$r <- c(.45,.10)
th_y_sim_homo$b <- c(1.2,-.5)
th_y_sim_homo$s <- c(.01,.02,.05,.04)

sim_homo <- simPowLawMixIndep(th_y_sim_homo,th_x_sim,100,modSpec_homo)
expect_error(
  fit_homo <- fit_theta_y(sim_homo$x,sim_homo$Y,modSpec_homo),
  NA
)

# Multi-variable, heteroskedastic (sd_x), NA in hetGroups
modSpec_heterox <- list(meanSpec='powLaw')
modSpec_heterox$J <- 2
modSpec_heterox$K <- 2
modSpec_heterox$M <- c(2,2)
modSpec_heterox$hetSpec <- 'sd_x'
modSpec_heterox$hetGroups <- c(NA,2,1,1)
modSpec_heterox$cdepSpec <- 'indep'

th_y_sim_heterox <- list(modSpec=modSpec_heterox)
th_y_sim_heterox$rho <- c(.75,.5)
th_y_sim_heterox$tau <- list()
th_y_sim_heterox$tau[[1]] <- c(2,5)
th_y_sim_heterox$tau[[2]] <- c(1,3)
th_y_sim_heterox$a <- c(2,3)
th_y_sim_heterox$r <- c(.45,.10)
th_y_sim_heterox$b <- c(1.2,-.5)
th_y_sim_heterox$s <- c(.01,.02,.05,.04)
th_y_sim_heterox$kappa <- kappa_full[1:2]


sim_heterox <- simPowLawMixIndep(th_y_sim_heterox,th_x_sim,100,modSpec_heterox)
expect_error(
  fit_heterox <- fit_theta_y(sim_heterox$x,sim_heterox$Y,modSpec_heterox),
  NA
)

# Multi-variable, heteroskedastic (sd_resp), NA in hetGroups
modSpec_heteror <- list(meanSpec='powLaw')
modSpec_heteror$J <- 2
modSpec_heteror$K <- 2
modSpec_heteror$M <- c(2,2)
modSpec_heteror$hetSpec <- 'sd_resp'
modSpec_heteror$hetGroups <- c(NA,2,1,1)
modSpec_heteror$cdepSpec <- 'indep'

th_y_sim_heteror <- list(modSpec=modSpec_heteror)
th_y_sim_heteror$rho <- c(.75,.5)
th_y_sim_heteror$tau <- list()
th_y_sim_heteror$tau[[1]] <- c(2,5)
th_y_sim_heteror$tau[[2]] <- c(1,3)
th_y_sim_heteror$a <- c(2,3)
th_y_sim_heteror$r <- c(.45,.10)
th_y_sim_heteror$b <- c(1.2,-.5)
th_y_sim_heteror$s <- c(.01,.02,.05,.04)
th_y_sim_heteror$kappa <- kappa_full[1:2]

sim_heteror <- simPowLawMixIndep(th_y_sim_heteror,th_x_sim,N,modSpec_heteror)
expect_error(
  fit_heteror <- fit_theta_y(sim_heteror$x,sim_heteror$Y,modSpec_heteror),
  NA
)

# Multi-variable, heteroskedastic (sd_x), single heteroskedastic parameter
modSpec_heterox <- list(meanSpec='powLaw')
modSpec_heterox$J <- 2
modSpec_heterox$K <- 2
modSpec_heterox$M <- c(2,2)
modSpec_heterox$hetSpec <- 'sd_x'
modSpec_heterox$hetGroups <- c(1,1,1,1)
modSpec_heterox$cdepSpec <- 'indep'

th_y_sim_heterox <- list(modSpec=modSpec_heterox)
th_y_sim_heterox$rho <- c(.75,.5)
th_y_sim_heterox$tau <- list()
th_y_sim_heterox$tau[[1]] <- c(2,5)
th_y_sim_heterox$tau[[2]] <- c(1,3)
th_y_sim_heterox$a <- c(2,3)
th_y_sim_heterox$r <- c(.45,.10)
th_y_sim_heterox$b <- c(1.2,-.5)
th_y_sim_heterox$s <- c(.01,.02,.05,.04)
th_y_sim_heterox$kappa <- kappa_full[1]


sim_heterox <- simPowLawMixIndep(th_y_sim_heterox,th_x_sim,100,modSpec_heterox)
expect_error(
  fit_heterox <- fit_theta_y(sim_heterox$x,sim_heterox$Y,modSpec_heterox),
  NA
)

# Multi-variable, heteroskedastic (sd_resp), single heteroskedastic parameter
modSpec_heteror <- list(meanSpec='powLaw')
modSpec_heteror$J <- 2
modSpec_heteror$K <- 2
modSpec_heteror$M <- c(2,2)
modSpec_heteror$hetSpec <- 'sd_resp'
modSpec_heteror$hetGroups <- c(1,1,1,1)
modSpec_heteror$cdepSpec <- 'indep'

th_y_sim_heteror <- list(modSpec=modSpec_heteror)
th_y_sim_heteror$rho <- c(.75,.5)
th_y_sim_heteror$tau <- list()
th_y_sim_heteror$tau[[1]] <- c(2,5)
th_y_sim_heteror$tau[[2]] <- c(1,3)
th_y_sim_heteror$a <- c(2,3)
th_y_sim_heteror$r <- c(.45,.10)
th_y_sim_heteror$b <- c(1.2,-.5)
th_y_sim_heteror$s <- c(.01,.02,.05,.04)
th_y_sim_heteror$kappa <- kappa_full[1]

sim_heteror <- simPowLawMixIndep(th_y_sim_heteror,th_x_sim,N,modSpec_heteror)
expect_error(
  fit_heteror <- fit_theta_y(sim_heteror$x,sim_heteror$Y,modSpec_heteror),
  NA
)

# Multi-variable, heteroskedastic (sd_x), multiple heteroskedastic parameters
modSpec_heterox <- list(meanSpec='powLaw')
modSpec_heterox$J <- 2
modSpec_heterox$K <- 2
modSpec_heterox$M <- c(2,2)
modSpec_heterox$hetSpec <- 'sd_x'
modSpec_heterox$hetGroups <- c(1,2,1,1)
modSpec_heterox$cdepSpec <- 'indep'

th_y_sim_heterox <- list(modSpec=modSpec_heterox)
th_y_sim_heterox$rho <- c(.75,.5)
th_y_sim_heterox$tau <- list()
th_y_sim_heterox$tau[[1]] <- c(2,5)
th_y_sim_heterox$tau[[2]] <- c(1,3)
th_y_sim_heterox$a <- c(2,3)
th_y_sim_heterox$r <- c(.45,.10)
th_y_sim_heterox$b <- c(1.2,-.5)
th_y_sim_heterox$s <- c(.01,.02,.05,.04)
th_y_sim_heterox$kappa <- kappa_full[1:3]


sim_heterox <- simPowLawMixIndep(th_y_sim_heterox,th_x_sim,100,modSpec_heterox)
expect_error(
  fit_heterox <- fit_theta_y(sim_heterox$x,sim_heterox$Y,modSpec_heterox),
  NA
)

# Multi-variable, heteroskedastic (sd_resp), multiple heteroskedastic parameters
modSpec_heteror <- list(meanSpec='powLaw')
modSpec_heteror$J <- 2
modSpec_heteror$K <- 2
modSpec_heteror$M <- c(2,2)
modSpec_heteror$hetSpec <- 'sd_resp'
modSpec_heteror$hetGroups <- c(1,2,1,1)
modSpec_heteror$cdepSpec <- 'indep'

th_y_sim_heteror <- list(modSpec=modSpec_heteror)
th_y_sim_heteror$rho <- c(.75,.5)
th_y_sim_heteror$tau <- list()
th_y_sim_heteror$tau[[1]] <- c(2,5)
th_y_sim_heteror$tau[[2]] <- c(1,3)
th_y_sim_heteror$a <- c(2,3)
th_y_sim_heteror$r <- c(.45,.10)
th_y_sim_heteror$b <- c(1.2,-.5)
th_y_sim_heteror$s <- c(.01,.02,.05,.04)
th_y_sim_heteror$kappa <- kappa_full[1:3]

sim_heteror <- simPowLawMixIndep(th_y_sim_heteror,th_x_sim,N,modSpec_heteror)
expect_error(
  fit_heteror <- fit_theta_y(sim_heteror$x,sim_heteror$Y,modSpec_heteror),
  NA
)

