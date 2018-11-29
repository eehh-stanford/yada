context("Siler")

# a0 is the baseline parameter vector used to test Siler functions. It is from
# Gage and Dyke 1986, Table 2, Level 15.
a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)

test_that("test Siler inverse", {
  qin <- seq(0, .999, by = .001) # quantiles in
  x <- qsiler(qin, a0) # x from quantiles
  qout <- psiler(x, a0) # quantiles from x
  expect_equal(qin, qout)
})

test_that("test Siler random sampler", {
  set.seed(825201) # from random.org
  N <- 1000000
  x <- seq(0, 120, by = .01)
  cdf <- psiler(x, a0)
  xsamp <- rsiler(N, a0)
  expect_equal(cdf, ecdf(xsamp)(x), tolerance = 1e-3)
})

x0 <- 10
test_that("test Siler inverse with offset", {
  qin <- seq(0, .999, by = .001) # quantiles in
  x <- qsiler(qin, a0, x0) # x from quantiles
  qout <- psiler(x, a0, x0) # quantiles from x
  expect_equal(qin, qout)
})

test_that("test Siler random sampler with offset", {
  set.seed(400532) # from random.org
  N <- 1000000
  x <- seq(x0, 120, by = .01)
  cdf <- psiler(x, a0, x0)
  xsamp <- rsiler(N, a0, x0)
  expect_equal(cdf, ecdf(xsamp)(x), tolerance = 1e-3)
})

test_that("test Siler maximum likelihood estimation", {
  set.seed(291542) # from random.org
  N <- 1000000
  xsamp <- rsiler(N, a0)
  a1 <- a0 * runif(5, min = .9, max = 1.1) # jitter start
  silerFit <- fitSilerMaxLik(xsamp, a1, calcHessian = TRUE)
  expect_equal(a0, silerFit$a, tol = 1e-3)
})
