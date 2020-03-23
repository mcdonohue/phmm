
skip_if(grepl('i386', sessionInfo()$platform), 
  message = 'Skip tests on i386 platforms.')

test_that("phmm with random intercept as expected", {
  # coxph(Surv(time, status) ~ age + frailty(inst, df=4), lung)
  set.seed(20200316)
  fit.phmm <- phmm(Surv(time, status) ~ age + (1|inst),
    data=survival::lung, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
    NINIT = 10, MAXSTEP = 100, CONVERG=90)
  expect_equivalent(fit.phmm$coefficients, 0.01877795, tolerance = .001)
  expect_equivalent(fit.phmm$Sigma, 0.004011077, tolerance = .001)
})

test_that("phmm with two random effects fit as expected", {
  set.seed(20200316)
  fit2.phmm <- phmm(Surv(time, status) ~ age + (sex|inst),
    data=survival::lung, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
    NINIT = 10, MAXSTEP = 100, CONVERG=90)
  expect_equivalent(fit2.phmm$coefficients, 0.01891319, tolerance = .001)
  expect_equivalent(diag(fit2.phmm$Sigma), c(0.004186401, 0.002143200), tolerance = .002)
})
