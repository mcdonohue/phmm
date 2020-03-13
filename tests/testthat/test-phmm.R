context("phmm")

n <- 50      # total sample size
nclust <- 5  # number of clusters
clusters <- rep(1:nclust,each=n/nclust)
beta0 <- c(1,2)
set.seed(13)
#generate phmm data set
Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
    Z2=sample(0:1,n,replace=TRUE),
    Z3=sample(0:1,n,replace=TRUE))
b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
Wb <- matrix(0,n,2)
for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
Wb <- apply(Wb,1,sum)
T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
C <- runif(n,0,1)
time <- ifelse(T<C,T,C)
event <- ifelse(T<=C,1,0)
mean(event)
phmmd <- data.frame(Z)
phmmd$cluster <- clusters
phmmd$time <- time
phmmd$event <- event

fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster),
    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
    NINIT = 10, MAXSTEP = 200, CONVERG=90)

test_that("phmm fit as expected", {
    expect_equal(fit.phmm$Sigma[1,1], 0.005107908, tolerance = .001)
    expect_equal(fit.phmm$Sigma[2,2], 1.057397, tolerance = .02)
    expect_equivalent(fit.phmm$coefficients[1], 1.6084942, tolerance = .005)
    expect_equivalent(fit.phmm$coefficients[2], 0.5784425, tolerance = .03)
})
