\name{ loglik.cond}
\alias{ loglik.cond}
\title{
PHMM conditional log-likelihood
}
\description{
Function for computing log-likelihood conditional on the estimated random effects from an object of class \code{phmm} returned by \code{phmm}.
}
\usage{
 loglik.cond(x)
}
\arguments{
\item{x}{an object of class \code{phmm}.}
}
\value{
The PHMM log-likelihood conditional on the estimated random effects.}
\seealso{
\code{\link{phmm}}, \code{\link{phmm.cond.loglik}}}
\keyword{survival}