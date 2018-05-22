context("Siler")

# a0 is the baseline parameter vector used to test Siler functions. It is from
# Gage and Dyke 1986, Table 2, Level 15 Table 2. The article does not describe
# things well, but it seems that the last three parameters must be modified by
# multiplying by .01, .001, and .1, respectively.
a0 <- c(.175,1.40,.368*.01,.075*.001,.917*.1) 

test_that("test Siler inverse", {
    qin  <- seq(0,.999,by=.001) # quantiles in
    x    <- qsiler(qin,a0) # x from quantiles
    qout <- psiler(x,a0) # quantiles from x
    expect_equal(qin,qout)
})

test_that("test Siler random sampler", {
    set.seed(825201) # from random.org
    N <- 1000000
    x <- seq(0,120,by=.01)
    cdf <- psiler(x,a0)
    xsamp <- rsiler(N,a0)
    expect_equal(cdf,ecdf(xsamp)(x),tolerance=1e-3)
})
