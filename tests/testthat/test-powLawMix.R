# Test the functions in powLawMix

# From random.org between 1 and 1,000,000
set.seed(937974)

# Simulate data with two ordinal and two continuous variables with no
# correlations
rho <- c(.75,.5)
tau <- list()
tau[[1]] <- c(2,5)
tau[[2]] <- c(1,3)
a <- c(2,3)
r <- c(.45,.10)
b <- c(1.2,-.5)
s <- c(.01,.02,.05,.04)
kappa_full <- c(.005,.0075,.01,.0125)

# Test a homoskedastic model with two ordinal and two continuous variables
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec = 'none'
modSpec$cdepSpec = 'indep'

th_y <- c(rho,unlist(tau),a,r,b,s)

# Test extract_th_v for both ordinal variables
expect_equal(
  extract_th_v(th_y,modSpec,1),
  c(rho[1],tau[[1]],s[1])
)

expect_equal(
  extract_th_v(th_y,modSpec,2),
  c(rho[2],tau[[2]],s[2])
)

# Test extract_th_w for both continuos variables
expect_equal(
  extract_th_w(th_y,modSpec,1),
  c(a[1],r[1],b[1],s[3])
)

expect_equal(
  extract_th_w(th_y,modSpec,2),
  c(a[2],r[2],b[2],s[4])
)

# Test a heteroskedastic model (sd_x) with two ordinal and two continuous variables
modSpec<- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec = 'sd_x'
modSpec$hetGroups <- 1:4
modSpec$cdepSpec = 'indep'

th_y <- c(rho,unlist(tau),a,r,b,s,kappa_full)

# Test extract_th_v for both ordinal variables
expect_equal(
  extract_th_v(th_y,modSpec,1),
  c(rho[1],tau[[1]],s[1],kappa_full[1])
)

expect_equal(
  extract_th_v(th_y,modSpec,2),
  c(rho[2],tau[[2]],s[2],kappa_full[2])
)

# Test extract_th_w for both continuos variables
expect_equal(
  extract_th_w(th_y,modSpec,1),
  c(a[1],r[1],b[1],s[3],kappa_full[3])
)

expect_equal(
  extract_th_w(th_y,modSpec,2),
  c(a[2],r[2],b[2],s[4],kappa_full[4])
)

# Test a heteroskedastic model (sd_resp) with two ordinal and two continuous variables
modSpec<- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec = 'sd_resp'
modSpec$hetGroups <- 1:4
modSpec$cdepSpec = 'indep'

th_y <- c(rho,unlist(tau),a,r,b,s,kappa_full)

# Test extract_th_v for both ordinal variables
expect_equal(
  extract_th_v(th_y,modSpec,1),
  c(rho[1],tau[[1]],s[1],kappa_full[1])
)

expect_equal(
  extract_th_v(th_y,modSpec,2),
  c(rho[2],tau[[2]],s[2],kappa_full[2])
)

# Test extract_th_w for both continuos variables
expect_equal(
  extract_th_w(th_y,modSpec,1),
  c(a[1],r[1],b[1],s[3],kappa_full[3])
)

expect_equal(
  extract_th_w(th_y,modSpec,2),
  c(a[2],r[2],b[2],s[4],kappa_full[4])
)

# Test a heteroskedastic model (sd_x) with two ordinal and two continuous variables
# and for which hetGroups is not 1:4
modSpec<- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec = 'sd_x'
modSpec$hetGroups <- c(2,NA,1,2)
modSpec$cdepSpec = 'indep'

th_y <- c(rho,unlist(tau),a,r,b,s,kappa_full[1],kappa_full[2])

# Test extract_th_v for both ordinal variables
expect_equal(
  extract_th_v(th_y,modSpec,1),
  c(rho[1],tau[[1]],s[1],kappa_full[2])
)

expect_equal(
  extract_th_v(th_y,modSpec,2),
  c(rho[2],tau[[2]],s[2],0)
)

# Test extract_th_w for both continuos variables
expect_equal(
  extract_th_w(th_y,modSpec,1),
  c(a[1],r[1],b[1],s[3],kappa_full[1])
)

expect_equal(
  extract_th_w(th_y,modSpec,2),
  c(a[2],r[2],b[2],s[4],kappa_full[2])
)

# Test a heteroskedastic model (sd_resp) with two ordinal and two continuous variables
# and for which hetGroups is not 1:4
modSpec<- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec = 'sd_resp'
modSpec$hetGroups <- c(2,NA,1,2)
modSpec$cdepSpec = 'indep'

th_y <- c(rho,unlist(tau),a,r,b,s,kappa_full[1],kappa_full[2])

# Test extract_th_v for both ordinal variables
expect_equal(
  extract_th_v(th_y,modSpec,1),
  c(rho[1],tau[[1]],s[1],kappa_full[2])
)

expect_equal(
  extract_th_v(th_y,modSpec,2),
  c(rho[2],tau[[2]],s[2],0)
)

# Test extract_th_w for both continuos variables
expect_equal(
  extract_th_w(th_y,modSpec,1),
  c(a[1],r[1],b[1],s[3],kappa_full[1])
)

expect_equal(
  extract_th_w(th_y,modSpec,2),
  c(a[2],r[2],b[2],s[4],kappa_full[2])
)

N <- 100
xmin <- 0
xmax <- 20
th_x_sim <- list(fitType='uniform',xmin=xmin,xmax=xmax)

# Check simulation for homoskedastic case
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

expect_error(
  sim_homo <- simPowLawMixIndep(th_y_sim_homo,th_x_sim,100,modSpec_homo),
  NA
)

expect_equal(
  names(sim_homo),
  c('x','Ystar','Y')
)

expect_equal(
  length(sim_homo$x),
  N
)

expect_equal(
  any(is.na(sim_homo$x)),
  F
)

expect_equal(
  dim(sim_homo$Ystar),
  c(4,N)
)

expect_equal(
  any(is.na(sim_homo$Ystar)),
  F
)

expect_equal(
  dim(sim_homo$Y),
  c(4,N)
)

expect_equal(
  any(is.na(sim_homo$Y)),
  F
)

# Check simulation for heteroskedastic case (sd_x)
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

expect_error(
  sim_heterox <- simPowLawMixIndep(th_y_sim_heterox,th_x_sim,N,modSpec_heterox),
  NA
)

expect_equal(
  names(sim_heterox),
  c('x','Ystar','Y')
)

expect_equal(
  length(sim_heterox$x),
  N
)

expect_equal(
  any(is.na(sim_heterox$x)),
  F
)

expect_equal(
  dim(sim_heterox$Ystar),
  c(4,N)
)

expect_equal(
  any(is.na(sim_heterox$Ystar)),
  F
)

expect_equal(
  dim(sim_heterox$Y),
  c(4,N)
)

expect_equal(
  any(is.na(sim_heterox$Y)),
  F
)

# Check simulation for heteroskedastic case (sd_resp)
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

expect_error(
  sim_heteror <- simPowLawMixIndep(th_y_sim_heteror,th_x_sim,N,modSpec_heteror),
  NA
)

expect_equal(
  names(sim_heteror),
  c('x','Ystar','Y')
)

expect_equal(
  length(sim_heteror$x),
  N
)

expect_equal(
  any(is.na(sim_heteror$x)),
  F
)

expect_equal(
  dim(sim_heteror$Ystar),
  c(4,N)
)

expect_equal(
  any(is.na(sim_heteror$Ystar)),
  F
)

expect_equal(
  dim(sim_heteror$Y),
  c(4,N)
)

expect_equal(
  any(is.na(sim_heteror$Y)),
  F
)

