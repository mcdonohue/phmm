context("phmm")

set.seed(20200316)
fit.phmm <- phmm(Surv(time, status) ~ age + (1|inst),
  data=survival::lung, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
  NINIT = 10, MAXSTEP = 100, CONVERG=90)

test_that("phmm with random intercept as expected", {
    expect_equivalent(fit.phmm$coefficients, 0.01877795, tolerance = .001)
})

# coxph(Surv(time, status) ~ age + frailty(inst, df=4), lung)