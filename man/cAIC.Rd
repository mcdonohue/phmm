\name{cAIC}
\alias{cAIC}
\alias{cAIC.coxph}
\alias{cAIC.phmm}
\title{Conditional Akaike Information Criterion for PHMM}
\usage{
  \method{cAIC}{phmm}(object, method = "direct", ..., k = 2)
  \method{cAIC}{coxph}(object, ..., k = 2)
}
\arguments{
  \item{object}{A fitted PHMM model object of class
  \code{phmm}.}

  \item{method}{Passed to \code{\link{traceHat}}. Options
  include "direct", "pseudoPois", or "HaLee". The methods
  "direct" and "HaLee" are algebraically equivalent.}

  \item{...}{Optionally more fitted model objects.}

  \item{k}{numeric, the penalty per parameter to be used;
  the default k = 2 conforms with the classical AIC.}
}
\value{
  Returns a numeric value of the cAIC corresonding to the
  PHMM fit.
}
\description{
  Function calculating the conditional Akaike information
  criterion (Vaida and Blanchard 2005) for PHMM fitted
  model objects, according to the formula
  \eqn{-2*log-likelihood + k*rho}{-2*log-likelihood +
  k*rho}, where \eqn{rho}{rho} represents the "effective
  degrees of freedom" in the sense of Hodges and Sargent
  (2001). The function uses the log-likelihood conditional
  on the estimated random effects; and trace of the "hat
  matrix", using the generalized linear mixed model
  formulation of PHMM, to estimate \eqn{rho}{rho}. The
  default k = 2, conforms with the usual AIC.
}
\references{
  Vaida, F, and Blanchard, S. 2005. Conditional Akaike
  information for mixed-effects models. Biometrika, 92(2),
  351-.

  Donohue, MC, Overholser, R, Xu, R, and Vaida, F (January
  01, 2011). Conditional Akaike information under
  generalized linear and proportional hazards mixed models.
  \emph{Biometrika}, 98, 3, 685-700.

  Breslow, NE, Clayton, DG. (1993). Approximate Inference
  in Generalized Linear Mixed Models. Journal of the
  American Statistical Association, Vol. 88, No. 421, pp.
  9-25.

  Whitehead, J. (1980). Fitting Cox\'s Regression Model to
  Survival Data using GLIM. Journal of the Royal
  Statistical Society. Series C, Applied statistics, 29(3),
  268-.

  Hodges, JS, and Sargent, DJ. 2001. Counting degrees of
  freedom in hierarchical and other richly-parameterised
  models. Biometrika, 88(2), 367-.
}
\seealso{
  \code{\link{phmm}}, \code{\link[stats]{AIC}}
}
\keyword{survival}
\examples{
\dontrun{
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

summary(fit.lmer)
}
}
