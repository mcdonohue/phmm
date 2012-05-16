\name{ phmm.cond.loglik}
\alias{ phmm.cond.loglik}
\title{
PHMM conditional log-likelihood
}
\description{
Function for computing log-likelihood conditional on the estimated random effects from the data and specified parameter estimates of a PHMM.
}
\usage{
 phmm.cond.loglik(time, delta, z, beta, w, b)
}
\arguments{
  \item{time}{Follow-up time (right censored data).}
  \item{delta}{The status indicator (0=alive, 1=dead; or \code{TRUE}=dead, \code{FALSE}=alive).}
  \item{z}{Numeric matrix (\code{N}x\code{nfixed}) of covariates for fixed effects.}
  \item{beta}{Fitted fixed effects coefficients (\code{p}-vector).}
  \item{w}{Numeric matrix (\code{N}x\code{nrandom}) of covariates for random effects.}
  \item{b}{Numeric matrix (\code{N}x\code{nrandom}) of random effects estimates.}
}
\value{
The PHMM log-likelihood conditional on the estimated random effects.}
\seealso{
\code{\link{phmm}}, \code{\link{loglik.cond}}}
\keyword{survival}