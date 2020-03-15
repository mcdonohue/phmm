context("phmm")

fit.phmm <- phmm(Surv(time, status) ~ age + (1|inst),
  data=survival::lung, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
  NINIT = 10, MAXSTEP = 100, CONVERG=90)

test_that("phmm with random intercept as expected", {
    expect_equivalent(fit.phmm$Sigma, 0.004059414, tolerance = .002)
    expect_equivalent(fit.phmm$coefficients, 0.0187804, tolerance = .01)
})
