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
z <- c(.25,-.2,.3,.5,.05,.15)
kappa1 <- .01
kappa2 <- .02

x1 <- 1
x2 <- 2
#modSpec <- list(meanSpec='powLaw')
#modSpec$J <- 2
#modSpec$K <- 2
#modSpec$M <- c(2,2)
#modSpec$hetSpec  <- 'none'
#modSpec$cdepSpec <- 'none'

## A single ordinal variable

# Test a model with only the one ordinal variable
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

# getIntegInfo_theta_y
expect_equal(
  getIntegInfo_theta_y(theta_y_list,0),
  list(doIntegral=T,limArray=t(matrix(c(-Inf,tau1[1]))))
)

expect_equal(
  getIntegInfo_theta_y(theta_y_list,1),
  list(doIntegral=T,limArray=t(matrix(c(tau1[1],tau1[2]))))
)

expect_equal(
  getIntegInfo_theta_y(theta_y_list,2),
  list(doIntegral=T,limArray=t(matrix(c(tau1[2],Inf))))
)

# calcLogLik_theta_y
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

expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,2),
  log(1-pnorm((tau1[2]-x1^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,2),
  log(1-pnorm((tau1[2]-x2^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),t(matrix(c(0,1)))),
  log(pnorm((tau1[1]-x1^rho1)/s1)) + log(pnorm((tau1[2]-x2^rho1)/s1) - pnorm((tau1[1]-x2^rho1)/s1))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),t(matrix(c(NA,2)))),
  log(1-pnorm((tau1[2]-x2^rho1)/s1))
)

# Test a model with only the one continuous variable
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 0
modSpec$K <- 1
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'none'

theta_y_list <- list(modSpec=modSpec)
theta_y_list$a <- a1
theta_y_list$r <- r1
theta_y_list$b <- b1
theta_y_list$s <- s3


# calcLogLik_theta_y
w1 <- 1.2
w2 <- 1.4
expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,w1),
  dnorm(w1,a1*x1^r1+b1,s3,log=T)
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,w2),
  dnorm(w2,a1*x2^r1+b1,s3,log=T)
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),t(matrix(c(w1,w2)))),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + dnorm(w2,a1*x2^r1+b1,s3,log=T)
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),t(matrix(c(NA,w2)))),
  dnorm(w2,a1*x2^r1+b1,s3,log=T)
)

# Test a model with one ordinal and one continuous variable
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 1
modSpec$M <- 2
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'dep'
modSpec$cdepGroups <- c(1,2)

theta_y_list <- list(modSpec=modSpec)
theta_y_list$rho <- rho1
theta_y_list$tau <- tau1
theta_y_list$a <- a1
theta_y_list$r <- r1
theta_y_list$b <- b1
theta_y_list$s <- c(s1,s3)
theta_y_list$z <- z[2]

w1 <- 1.2
w2 <- 1.4
# getIntegInfo_theta_y
expect_equal(
  getIntegInfo_theta_y(theta_y_list,c(0,w1)),
  list(doIntegral=c(T,F),limArray=matrix(c(-Inf,NA,tau1[1],NA),nrow=2))
)

expect_equal(
  getIntegInfo_theta_y(theta_y_list,c(NA,w1)),
  list(doIntegral=c(T,F),limArray=matrix(c(-Inf,NA,Inf,NA),nrow=2))
)

expect_equal(
  getIntegInfo_theta_y(theta_y_list,c(NA,NA)),
  list(doIntegral=c(T,T),limArray=matrix(c(-Inf,-Inf,Inf,Inf),nrow=2))
)

# calcLogLik_theta_y
mu_bar1 <- x1^rho1 + z[2]*s1/s3*(w1-a1*x1^r1-b1)
s1_bar <- s1*sqrt(1-z[2]^2)
mu_bar2 <- x2^rho1 + z[2]*s1/s3*(w2-a1*x2^r1-b1)
expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,matrix(c(0,w1))),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + log(pnorm((tau1[1]-mu_bar1)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,matrix(c(1,w1))),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + log(pnorm((tau1[2]-mu_bar1)/s1_bar) - pnorm((tau1[1]-mu_bar1)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,matrix(c(2,w1))),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + log(1 - pnorm((tau1[2]-mu_bar1)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,matrix(c(0,w2))),
  dnorm(w2,a1*x2^r1+b1,s3,log=T) + log(pnorm((tau1[1]-mu_bar2)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,matrix(c(1,w2))),
  dnorm(w2,a1*x2^r1+b1,s3,log=T) + log(pnorm((tau1[2]-mu_bar2)/s1_bar) - pnorm((tau1[1]-mu_bar2)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,matrix(c(2,w2))),
  dnorm(w2,a1*x2^r1+b1,s3,log=T) + log(1 - pnorm((tau1[2]-mu_bar2)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),matrix(c(0,w1,1,w2),nrow=2)),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + log(pnorm((tau1[1]-mu_bar1)/s1_bar)) + dnorm(w2,a1*x2^r1+b1,s3,log=T) + log(pnorm((tau1[2]-mu_bar2)/s1_bar) - pnorm((tau1[1]-mu_bar2)/s1_bar))
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),matrix(c(NA,w1,1,w2),nrow=2)),
  dnorm(w1,a1*x1^r1+b1,s3,log=T) + dnorm(w2,a1*x2^r1+b1,s3,log=T) + log(pnorm((tau1[2]-mu_bar2)/s1_bar) - pnorm((tau1[1]-mu_bar2)/s1_bar))
)

# Test a model with two ordinal and two continuous variable
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'dep'
modSpec$cdepGroups <- c(1,2,3,4)

theta_y_list <- list(modSpec=modSpec)
theta_y_list$rho <- c(rho1,rho2)
theta_y_list$tau <- list()
theta_y_list$tau[[1]] <- tau1
theta_y_list$tau[[2]] <- tau2
theta_y_list$a <- c(a1,a2)
theta_y_list$r <- c(r1,r2)
theta_y_list$b <- c(b1,b2)
theta_y_list$s <- c(s1,s2,s3,s4)
theta_y_list$z <- z

Y <- matrix(c(0,1,100,50,1,2,160,92),nrow=4)
# Directly Calculate likelihood for both observations
# Ordinal means:
g1 <- x1^c(rho1,rho2)
g2 <- x2^c(rho1,rho2)
h1 <- c(a1,a2)*x1^c(r1,r2) + c(b1,b2)
h2 <- c(a1,a2)*x2^c(r1,r2) + c(b1,b2)
S1 <- get_Sigma(theta_y_list,x1)
S2 <- get_Sigma(theta_y_list,x2)

condMean1 <- g1 + S1[1:2,3:4] %*% solve(S1[3:4,3:4]) %*% as.matrix(c(100,50) - h1)
condCov1 <- S1[1:2,1:2] - S1[1:2,3:4] %*% solve(S1[3:4,3:4]) %*% S1[3:4,1:2]
lik1 <- mvtnorm::dmvnorm(c(100,50),h1,S1[3:4,3:4],log=T)
condIntegral1 <- mvtnorm::pmvnorm(lower=c(-Inf,tau2[1]),upper=c(tau1[1],tau2[2]),mean=as.vector(condMean1),sigma=condCov1)
lik1 <- lik1 + log(as.numeric(condIntegral1))
expect_equal(
  calcLogLik_theta_y(theta_y_list,x1,Y[,1]),
  lik1
)

condMean2 <- g2 + S2[1:2,3:4] %*% solve(S2[3:4,3:4]) %*% as.matrix(c(160,92) - h2)
condCov2 <- S2[1:2,1:2] - S2[1:2,3:4] %*% solve(S2[3:4,3:4]) %*% S2[3:4,1:2]
lik2 <- mvtnorm::dmvnorm(c(160,92),h2,S2[3:4,3:4],log=T)
condIntegral2 <- mvtnorm::pmvnorm(lower=c(tau1[1],tau2[2]),upper=c(tau1[2],Inf),mean=as.vector(condMean2),sigma=condCov2)
lik2 <- lik2 + log(as.numeric(condIntegral2))
expect_equal(
  calcLogLik_theta_y(theta_y_list,x2,Y[,2]),
  lik2
)

expect_equal(
  calcLogLik_theta_y(theta_y_list,c(x1,x2),Y),
  lik1 + lik2
)

# For powLawMixIndepNegLogLik, test a conditionally independent model with two ordinal and two continuous variable
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,2)
modSpec$hetSpec  <- 'none'
modSpec$cdepSpec <- 'indep'

theta_y_list <- list(modSpec=modSpec)
theta_y_list$rho <- c(rho1,rho2)
theta_y_list$tau <- list()
theta_y_list$tau[[1]] <- tau1
theta_y_list$tau[[2]] <- tau2
theta_y_list$a <- c(a1,a2)
theta_y_list$r <- c(r1,r2)
theta_y_list$b <- c(b1,b2)
theta_y_list$s <- c(s1,s2,s3,s4)

Y <- matrix(c(0,1,100,50,1,2,160,92),nrow=4)
# Directly Calculate likelihood for both observations
# Ordinal means:
g1 <- x1^c(rho1,rho2)
g2 <- x2^c(rho1,rho2)
h1 <- c(a1,a2)*x1^c(r1,r2) + c(b1,b2)
h2 <- c(a1,a2)*x2^c(r1,r2) + c(b1,b2)

lik1_ord <- c(pnorm((tau1[1]-g1[1])/s1),pnorm((tau2[2]-g1[2])/s2) - pnorm((tau2[1]-g1[2])/s2))
lik1_cont <- dnorm(Y[3:4,1],h1,c(s3,s4))

lik2_ord <- c(pnorm((tau1[2]-g2[1])/s1) - pnorm((tau1[1]-g2[1])/s1),1 - pnorm((tau2[2]-g2[2])/s2))
lik2_cont <- dnorm(Y[3:4,2],h2,c(s3,s4))

eta <- -sum(log(c(lik1_ord,lik1_cont,lik2_ord,lik2_cont)))

theta_y_vect <- theta_y_list2vect(theta_y_list)
expect_equal(
  powLawMixIndepNegLogLik(theta_y_vect,c(x1,x2),Y,modSpec),
  eta
)

