## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
library(yada)
Dplot_MSR <- seq(-0.5, 0.5, len = 1001) # The locations at which to plot the PDFs
plot(Dplot_MSR, calcMSRDensityConstHaz(Dplot_MSR, .25, component = "event"), type = "l", xlim = c(-0.5, 0.5), ylim = c(0, 2), xlab = "Residual Value", ylab = "MSR PDF", lwd = 3, col = "gray", lty = 2)
lines(Dplot_MSR, calcMSRDensityConstHaz(Dplot_MSR, .25), lwd = 3, col = "gray")
lines(Dplot_MSR, calcMSRDensityConstHaz(Dplot_MSR, 0), lwd = 3, col = "black")
legend("topright", c("No Censoring", "25% Censoring"), lwd = 3, col = c("black", "gray"))

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
Dplot_NMR <- seq(-1, 2, len = 1001) # The locations at which to plot the PDFs
plot(Dplot_NMR, calcNMRDensityConstHaz(Dplot_NMR, .25, component = "event"), type = "l", xlim = c(-1, 2), ylim = c(0, 1), xlab = "Residual Value", ylab = "NMR PDF", lwd = 3, col = "gray", lty = 2)
lines(Dplot_NMR, calcNMRDensityConstHaz(Dplot_NMR, .25), lwd = 3, col = "gray")
lines(Dplot_NMR, calcNMRDensityConstHaz(Dplot_NMR, 0), lwd = 3, col = "black")
legend("topright", c("No Censoring", "25% Censoring"), lwd = 3, col = c("black", "gray"))

## ------------------------------------------------------------------------
# Set the random number seed
set.seed(60498) # Chosen to effectively illustrate the Kolmogorov-Smirnov statistic

n <- 4
rho <- .2
b1 <- 1 / 20
b2 <- b1 * rho / (1 - rho) # Censoring rate [rho is only apporoximately the censoring ratio given a linear event hazard]

a1 <- b1 / 10
survObjLin <- simLinHaz(n, b1, a1, b2)
print(survObjLin)

## ------------------------------------------------------------------------
fitConst <- fitConstHaz(survObjLin) # The wrong model
cumHazEvent <- fitConst$b1 * survObjLin[, "time"]
cumHazCens <- fitConst$b2 * survObjLin[, "time"]
cumHazTot <- cumHazEvent + cumHazCens
print(cumHazTot)

## ------------------------------------------------------------------------
# Martingale Survival Residual treating censoring as a competing risk
msrTot <- calcMSR(cumHazTot)
ksLin <- ks.test(msrTot, punif, -.5, .5) # The Kolmogorov-Smirnov test

## ------------------------------------------------------------------------
cdf <- calcStoppingCDF(msrTot, rep(F, length(msrTot)), where = "both")
indMax <- which.max(c(abs(cdf$y + .5 - cdf$before), abs(cdf$y + .5 - cdf$after)))
isBefore <- indMax <= length(cdf$y)
ind0 <- isBefore * (indMax) + (1 - isBefore) * (indMax - 4)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
library(shape)
plot(x = NULL, y = NULL, frame = T, xlim = c(-.5, .5), ylim = c(0, 1), lwd = 3, xlab = "Residual Value", ylab = "Empirical CDF")
lines(c(-.5, .5), c(0, 1), col = "gray", lwd = 3)
lines(c(1, 1) * cdf$y[ind0], c(0, cdf$before[ind0]), col = "gray", lwd = 3, lty = "dashed")
if (isBefore) {
  stop("This should not happen for the seed used")
} else {
  lines(c((cdf$y[ind0] - .5) / 2, cdf$y[ind0]), c(1, 1) * cdf$y[ind0] + .5, col = "gray", lwd = 3, lty = "dashed")
  lines(c((cdf$y[ind0] - .5) / 2, cdf$y[ind0]), c(1, 1) * cdf$after[ind0], col = "gray", lwd = 3, lty = "dashed")
  Arrows((cdf$y[ind0] - .5) / 2, cdf$y[ind0] + .5, (cdf$y[ind0] - .5) / 2, cdf$after[ind0], code = 3, col = "gray", arr.type = "triangle", arr.adj = -.5, arr.length = .25, arr.width = .25, lwd = 2)

  text(x = (cdf$y[ind0] - .5) / 2 - .1, y = (cdf$y[ind0] + .5 + cdf$after[ind0]) / 2, labels = expression(K[n]))
}
lines(c(-.5, cdf$y, .5), c(0, cdf$before, 1), type = "S", lwd = 3)

## ------------------------------------------------------------------------
print(cdf$y[ind0])

## ------------------------------------------------------------------------
print(cdf$after[ind0] - (cdf$y[ind0] + .5))

## ------------------------------------------------------------------------
ksTest <- ks.test(cdf$y, punif, -.5, .5)
print(ksTest$statistic)

## ------------------------------------------------------------------------
print(ksTest)

## ------------------------------------------------------------------------
# Constant transition intensities for event and censoring, with a correct
# candidate model that also has constant transition intensities.
set.seed(387341) # From random.org between 1 and 1,000,000

n <- 1000
rho <- .2
b1 <- 1 / 20 # Event rate
b2 <- b1 * rho / (1 - rho) # Censoring rate

# Create the simulated data
survObjConst <- simConstHaz(n, b1, b2)

## ------------------------------------------------------------------------
fitConst <- fitConstHaz(survObjConst)
cumHazEvent <- fitConst$b1 * survObjConst[, "time"]
cumHazCens <- fitConst$b2 * survObjConst[, "time"]
cumHazTot <- cumHazEvent + cumHazCens

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
msrCr <- calcMSR(cumHazTot)
hist(msrCr, 40, xlim = c(-1 / 2, 1 / 2), ylim = c(0, 2), freq = FALSE, xlab = "Residual Value", ylab = "Probability Density", main = NULL)
y <- seq(-1 / 2, 1 / 2, len = 1001)
lines(y, calcMSRDensityConstHaz(y, 0), lwd = 3)

## ------------------------------------------------------------------------
ksConst <- ks.test(msrCr, punif, -.5, .5) # The Kolmogorov-Smirnov test
print(ksConst$statistic)
print(ksConst$p.value)

## ------------------------------------------------------------------------
# Time dependent (linear) event intensity and constant censoring, with an incorrect
# candidate model that has constant transition intensities.
set.seed(51600) # From random.org between 1 and 1,000,000

n <- 1000
rho <- .2
b1 <- 1 / 20 # Event rate
b2 <- b1 * rho / (1 - rho) # Censoring rate [rho is only apporoximately the censoring ratio given a linear event hazard]
a1 <- b1 / 10

survObjLin <- simLinHaz(n, b1, a1, b2)

## ------------------------------------------------------------------------
fitConst <- fitConstHaz(survObjLin) # The wrong model
cumHazEvent <- fitConst$b1 * survObjLin[, "time"]
cumHazCens <- fitConst$b2 * survObjLin[, "time"]
cumHazTot <- cumHazEvent + cumHazCens

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
msrCr <- calcMSR(cumHazTot)
hist(msrCr, 40, xlim = c(-1 / 2, 1 / 2), ylim = c(0, 2), freq = FALSE, xlab = "Residual Value", ylab = "Probability Density", main = NULL)
y <- seq(-1 / 2, 1 / 2, len = 1001)
lines(y, calcMSRDensityConstHaz(y, 0), lwd = 3)

## ------------------------------------------------------------------------
ksLin <- ks.test(msrCr, punif, -.5, .5) # The Kolmogorov-Smirnov test
print(ksLin$statistic)
print(ksLin$p.value)

## ------------------------------------------------------------------------
alpha <- .05
nVect <- seq(10, 500, by = 10) # The vector of number of observations
powerVect <- rep(NA, length(nVect)) # The vector in which to store the power (1-beta)

# Nmc <- 50000
Nmc <- 500

options(warn = -1) # The numerical simulation yields very occasional ties. Suppress warnings this creates in ks.test
for (ii in 1:length(nVect)) {
  numSig <- 0
  n <- nVect[ii]
  for (m in 1:Nmc) {
    survObjLin <- simLinHaz(n, b1, a1, b2)
    fitConst <- fitConstHaz(survObjLin) # The wrong model
    cumHazEvent <- fitConst$b1 * survObjLin[, "time"]
    cumHazCens <- fitConst$b2 * survObjLin[, "time"]
    cumHazTot <- cumHazEvent + cumHazCens
    # Martingale Survival Residual treating censoring as a competing risk
    msrCr <- calcMSR(cumHazTot)
    ks2Pc <- ks.test(msrCr, punif, -.5, .5) # The Kolmogorov-Smirnov test
    if (ks2Pc$p.value <= alpha) {
      numSig <- numSig + 1
    }
    powerVect[ii] <- numSig / Nmc
  }
}
options(warn = 0) # Turn warnings back on

# Make the power zero at the origin for plotting
nVect <- c(0, nVect)
powerVect <- c(0, powerVect)

# Interpolate the curve's value for beta=.2 and beta = .1
y1 <- .8
y2 <- .9

x1 <- approx(powerVect, nVect, y1)$y
x2 <- approx(powerVect, nVect, y2)$y

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
plot(c(0, x1), c(y1, y1), xlab = "Sample Size", ylab = "Power", xlim = c(0, 500), ylim = c(0, 1), type = "l", col = "gray", lwd = 3, lty = "dashed")
lines(c(x1, x1), c(0, y1), col = "gray", lwd = 3, lty = "dashed")
lines(c(0, x2), c(y2, y2), col = "gray", lwd = 3, lty = "dashed")
lines(c(x2, x2), c(0, y2), col = "gray", lwd = 3, lty = "dashed")
lines(nVect, powerVect, lwd = 3)

## ------------------------------------------------------------------------
library(survival)

## Additive covariates incorrectly fit as proportional covariates
set.seed(2105984)

# Define the starting age (15 years)
x0 <- 15

# From Gage and Dyke 1986, Table 2, Level 15
afem <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1) # Female parameterization

# Males are assumed to have a higher hazard by an additive (constant) amount
maleBoost <- afem[3] * 2
amal <- afem + c(0, 0, maleBoost, 0, 0)

N <- 2000 # Number of simulated observations
# Sample randomly for male/female
female <- sample(c(TRUE, FALSE), N, replace = T)

# Sample from the relevant Siler hazard
x <- rep(NA, N) # vector of simulated ages
x[ female] <- rsiler(sum(female), afem, x0)
x[!female] <- rsiler(sum(!female), amal, x0)

## ------------------------------------------------------------------------
# Objective function for additive fit
fadd <- function(param, x, female, x0, a0) {
  # Calculate the objective function (negative log-likelihood)
  # given the input parameters (param), age-at-death vector (x),
  # and boolean vector indicating sex (female)
  a <- a0 * exp(param)
  afem <- a[1:5]
  amal <- afem + c(0, 0, a[6], 0, 0)

  xFemale <- x[female]
  xMale <- x[!female]

  nll <- nllsiler(afem, xFemale, x0) + nllsiler(amal, xMale, x0)
  return(nll)
}

# Objective function for multiplicative fit
fmult <- function(param, x, female, x0, a0) {
  # Calculate the objective function (negative log-likelihood)
  # given the input parameters (param), age-at-death vector (x),
  # and boolean vector indicating sex (female)
  a <- a0 * exp(param)
  afem <- a[1:5]
  amal <- afem * c(a[6], 1, a[6], 1, a[6])

  xFemale <- x[female]
  xMale <- x[!female]

  nll <- nllsiler(afem, xFemale, x0) + nllsiler(amal, xMale, x0)
  return(nll)
}

## ------------------------------------------------------------------------
# Fit the additive and multiplicative models and extract the fit parameters
a0_add <- c(afem, maleBoost)
paramFitAdd <- optim(rep(0, 6), fadd, x = x, female = female, x0 = x0, a0 = a0_add)
afem_add <- a0_add[1:5] * exp(paramFitAdd$par[1:5])
amal_add <- afem_add + c(0, 0, a0_add[6] * exp(paramFitAdd$par[6]), 0, 0)

a0_mult <- c(afem, 1)
paramFitMult <- optim(rep(0, 6), fmult, x = x, female = female, x0 = x0, a0 = a0_mult)
afem_mult <- a0_mult[1:5] * exp(paramFitMult$par[1:5])
amal_mult <- afem_mult * c(a0_mult[6] * exp(paramFitMult$par[6]), 1, a0_mult[6] * exp(paramFitMult$par[6]), 1, a0_mult[6] * exp(paramFitMult$par[6]))

## ------------------------------------------------------------------------
# Calculate the survivals, Martingale Survival Residuals (MSRs), and do the Kolmogorov-Smirnov tests
S_add <- rep(NA, N)
S_add[ female] <- ssiler(x[ female], afem_add, x0)
S_add[!female] <- ssiler(x[!female], amal_add, x0)

S_mult <- rep(NA, N)
S_mult[ female] <- ssiler(x[ female], afem_mult, x0)
S_mult[!female] <- ssiler(x[!female], amal_mult, x0)

Dhat_add <- 0.5 * (1 - 2 * S_add)
ks_add <- ks.test(Dhat_add, punif, -.5, .5) # The Kolmogorov-Smirnov test
ks_add_fem <- ks.test(Dhat_add[ female], punif, -.5, .5) # The Kolmogorov-Smirnov test
ks_add_mal <- ks.test(Dhat_add[!female], punif, -.5, .5) # The Kolmogorov-Smirnov test

Dhat_mult <- 0.5 * (1 - 2 * S_mult)
ks_mult <- ks.test(Dhat_mult, punif, -.5, .5) # The Kolmogorov-Smirnov test
ks_mult_fem <- ks.test(Dhat_mult[ female], punif, -.5, .5) # The Kolmogorov-Smirnov test
ks_mult_mal <- ks.test(Dhat_mult[!female], punif, -.5, .5) # The Kolmogorov-Smirnov test

## ------------------------------------------------------------------------
# Fit a Cox model, calculate the survivals / MSRs, and to the Kolmogorov-Smirnov tests
sampDf <- data.frame(x = x, female = female)
coxFit <- coxph(Surv(x) ~ female, data = sampDf)
cumHazCox <- as.numeric(predict(coxFit, sampDf, type = "expected"))
survCox <- exp(-cumHazCox)
Dhat_cox <- 0.5 * (1 - 2 * survCox)
ks_cox <- ks.test(Dhat_cox, punif, -.5, .5)
ks_cox_fem <- ks.test(Dhat_cox[ female], punif, -.5, .5) # The Kolmogorov-Smirnov test
ks_cox_mal <- ks.test(Dhat_cox[!female], punif, -.5, .5) # The Kolmogorov-Smirnov test

## ------------------------------------------------------------------------
# Create and print out a summary data frame
testName <- c("Add All", "Mlt All", "Cox All", "Add Fem", "Add Mal", "Mlt Fem", "Mlt Mal", "Cox Fem", "Cox Mal")
pvalVect <- c(ks_add$p.value, ks_mult$p.value, ks_cox$p.value, ks_add_fem$p.value, ks_add_mal$p.value, ks_mult_fem$p.value, ks_mult_mal$p.value, ks_cox_fem$p.value, ks_cox_mal$p.value)
ksStatVect <- c(ks_add$statistic, ks_mult$statistic, ks_cox$statistic, ks_add_fem$statistic, ks_add_mal$statistic, ks_mult_fem$statistic, ks_mult_mal$statistic, ks_cox_fem$statistic, ks_cox_mal$statistic)
dfOut <- data.frame(test = testName, statistic = ksStatVect, pval = pvalVect)
print(dfOut)

dfOut <- data.frame(test = testName, pval = pvalVect)
print(dfOut)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
# Plot the MSR histograms stratified by female/male
femCol <- rgb(1, 0, 0, .5)
malCol <- rgb(0, 1, 0, .5)
hist(Dhat_add[ female], 100, xlim = c(-1 / 2, 1 / 2), ylim = c(0, 2), freq = F, xlab = "Residual Value", ylab = "Probability Density", main = NULL, col = femCol)
hist(Dhat_add[!female], 100, add = T, freq = F, col = malCol)
y <- seq(-1 / 2, 1 / 2, len = 1001)
lines(y, calcMSRDensityConstHaz(y, 0), lwd = 3)
legend("topright", c("Female", "Male"), fill = c("red", "green"))
box()

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
hist(Dhat_mult[ female], 100, xlim = c(-1 / 2, 1 / 2), ylim = c(0, 2), freq = F, xlab = "Residual Value", ylab = "Probability Density", main = NULL, col = femCol)
hist(Dhat_mult[!female], 100, add = T, freq = F, col = malCol)
lines(y, calcMSRDensityConstHaz(y, 0), lwd = 3)
legend("topright", c("Female", "Male"), fill = c("red", "green"))
box()

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
hist(Dhat_cox[ female], 100, xlim = c(-1 / 2, 1 / 2), ylim = c(0, 2), freq = F, xlab = "Residual Value", ylab = "Probability Density", main = NULL, col = femCol)
hist(Dhat_cox[!female], 100, add = T, freq = F, col = malCol)
lines(y, calcMSRDensityConstHaz(y, 0), lwd = 3)
legend("topright", c("Female", "Male"), fill = c("red", "green"))
box()

