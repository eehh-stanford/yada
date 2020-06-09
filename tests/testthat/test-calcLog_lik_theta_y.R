# Test calcLogLik_theta_y
library(yada)
library(doParallel)
registerDoParallel(detectCores()-2)

rho1 <- .65
rho2 <- .75
tau1 <- c(1,1.5)
tau2 <- c(2,3)
a1   <- 110
a2   <- 40
r1   <- .45
r2   <- .55
b1   <- -45
b2   <- 15
s1 <- .25
s2 <- .5
s3 <- 20
s4 <- 10
z <- c(.25,-.2,.3,.5,.05)
kappa1 <- .01
kappa2 <- .02

x1 <- 1
x2 <- 2

# A single ordinal variable
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 0
modSpec$M <- 2
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'none'

theta_y_list <- list(modSpec=modSpec)
theta_y_list$rho <- rho1
theta_y_list$tau <- tau1
theta_y_list$s   <-   s1


modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'none'

# Test calculation with a scalar x for x1 and x2 and all three ordinal categories.
expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,0),
  log(pnorm((tau1[1]-x1^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,0),
  log(pnorm((tau1[1]-x2^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,1),
  log(pnorm((tau1[2]-x1^rho1)/s1) - pnorm((tau1[1]-x1^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,1),
  log(pnorm((tau1[2]-x2^rho1)/s1) - pnorm((tau1[1]-x2^rho1)/s1))
)
