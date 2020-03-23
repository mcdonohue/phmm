
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

test_that("phmm with three random effects fit as expected", {
  n <- 50      # total sample size
  nclust <- 5  # number of clusters
  clusters <- rep(1:nclust,each=n/nclust)
  beta0 <- c(1,2)
  set.seed(13)
  
  Z <-cbind(Z1=sample(0:1,n,replace=TRUE),
    Z2=sample(0:1,n,replace=TRUE),
    Z3=sample(0:1,n,replace=TRUE))
  b <- cbind(rep(rnorm(nclust), each=n/nclust),
    rep(rnorm(nclust), each=n/nclust))
  Wb <- matrix(0,n,2)
  for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
  Wb <- apply(Wb,1,sum)
  T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
  C <- runif(n,0,1)
  time <- ifelse(T<C,T,C)
  event <- ifelse(T <= C,1,0)
  phmmd <- data.frame(Z)
  phmmd$cluster <- clusters
  phmmd$time <- time
  phmmd$event <- event
  
  set.seed(20200316)
  fit2.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (Z1 + Z2|cluster), 
    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
    NINIT = 10, MAXSTEP = 100, CONVERG=90)
  expect_equivalent(fit2.phmm$coefficients, c(1.6169634, 0.5818235), tolerance = .001)
  # expect_equivalent(diag(fit2.phmm$Sigma), c(0.010257528, 0.006868375, 1.056108716), tolerance = .001)
})
