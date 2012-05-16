\name{ linear.predictors}
\alias{ linear.predictors}
\title{
PHMM Design
}
\description{
Internal function for extracting the linear predictors of the PHMM model from an object of class \code{phmm} returned by \code{phmm}.
}
\usage{
 linear.predictors(x)
}
\arguments{
\item{x}{an object of class \code{phmm}.}
}
\value{
A vector of estimates from call to \code{phmm} of \eqn{\beta'x_{ij}+w_{ij}'b_{i}}{beta'x[ij]+w[ij]'b[i]}.}
\seealso{
\code{\link{phmm}}}
\keyword{survival}