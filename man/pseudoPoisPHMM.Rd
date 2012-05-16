\name{ pseudoPoisPHMM}
\alias{ pseudoPoisPHMM}
\title{
Pseudo poisson data for fitting PHMM via GLMM
}
\description{
Function for generating a pseudo Poisson data set which can be used to fit a PHMM using GLMM software. This follows the mixed-model extension Whitehead (1980), who described how to fit Cox (fixed effects) models with GLM software.
}
\usage{
 pseudoPoisPHMM(x)
}
\arguments{
\item{x}{an object of class \code{phmm}.}
}
\value{
A \code{data.frame} with columns:
\item{time}{the event time;}
\item{N}{the number at risk at time \code{time};}
\item{m}{the number at risk (in the same cluster with same covariates) at time \code{time};}
\item{cluster}{the integer cluster indicator;}
\item{N}{the number at risk at time \code{time};}
\item{fixed effects covariates}{denoted \code{z1}, \code{z2}, etc.;}
\item{random effects covariates}{denoted \code{w1}, \code{w2}, etc.;}
\item{linear.predictors}{the linear predictors from the \code{phmm} fit (excluding the cumulative hazard estimates.}
}

\references{
Whitehead, J. (1980). Fitting Cox's Regression Model to Survival Data using GLIM. Journal of the Royal Statistical Society. Series C, Applied statistics, 29(3). 268-.}

\seealso{\code{\link{phmm}}, \code{\link{traceHat}}}
\examples{\dontrun{
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
T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]\%*\%beta0-Wb)
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
   NINIT = 10, MAXSTEP = 100, CONVERG=90)

# Same data can be fit with lmer,
# though the correlation structures are different.
poisphmmd <- pseudoPoisPHMM(fit.phmm)

library(lme4)
fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
  (-1+w1+w2|cluster)+offset(log(N)), 
  as.data.frame(as(poisphmmd, "matrix")), family=poisson)

fixef(fit.lmer)[c("z1","z2")]
fit.phmm$coef

VarCorr(fit.lmer)$cluster
fit.phmm$Sigma

logLik(fit.lmer)
fit.phmm$loglik

traceHat(fit.phmm)
}}
\keyword{survival}