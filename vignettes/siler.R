## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
a0 <- c(0.175, 1.40, 0.368 * 0.01, 0.075 * 0.001, 0.917 * 0.1)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
library(yada)
x <- seq(0, 100, by = .1)
hazSiler <- hsiler(x, a0)
th <- 3 # Line width for plots
plot(x, hazSiler, type = "l", xlab = "Age [years]", ylab = "Hazard", lwd = th)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
cumHazSiler <- chsiler(x, a0)
plot(x, cumHazSiler, type = "l", xlab = "Age [years]", ylab = "Cumulative Hazard", lwd = th)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
survSiler <- ssiler(x, a0)
plot(x, survSiler, type = "l", xlab = "Age [years]", ylab = "Survival", lwd = th)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
cdfSiler <- psiler(x, a0)
plot(x, cdfSiler, type = "l", xlab = "Age [years]", ylab = "Cumulative Density", lwd = th)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
pdfSiler <- dsiler(x, a0)
plot(x, pdfSiler, type = "l", xlab = "Age [years]", ylab = "Density", lwd = th)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
x0 <- 15
x2 <- seq(x0, 100, by = .1)
survSilerCont <- ssiler(x2, a0, x0)
plot(x2, survSilerCont, type = "l", xlab = "Age [years]", ylab = "Survival", lwd = th)
pdfSilerCont <- dsiler(x2, a0, x0)
plot(x2, pdfSilerCont, type = "l", xlab = "Age [years]", ylab = "Density", lwd = th)

## ------------------------------------------------------------------------
set.seed(896354) # from random.org
N <- 10000
xsamp <- rsiler(N, a0)

## ------------------------------------------------------------------------
a1 <- a0 * runif(5, min = .9, max = 1.1)
silerFit <- fitSilerMaxLik(xsamp, a1, calcHessian = TRUE)
print(silerFit$a)

## ------------------------------------------------------------------------
print(silerFit$fit)

## ----fig.align='center',fig.width=6,fig.height=6-------------------------
hist(xsamp, 100, freq = F, main = NA, xlab = "Age [years]", ylab = "Density")
lines(x, dsiler(x, a0), lty = 1, col = "1", lwd = th)
lines(x, dsiler(x, a1), lty = 2, col = "red", lwd = th)
lines(x, dsiler(x, silerFit$a), col = "red", lwd = th)
legend("topright", c("Target", "Initial", "Max Lik"), col = c("black", "red", "red"), lty = c(1, 2, 1))

## ------------------------------------------------------------------------
covMat <- solve(-silerFit$hessian) # invert the matrix

## ------------------------------------------------------------------------
standErr <- sqrt(diag(covMat))

## ------------------------------------------------------------------------
z <- silerFit$a / standErr

## ------------------------------------------------------------------------
pval <- pnorm(-abs(z))
print(pval)

