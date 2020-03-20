

#' Proportional Hazards with Mixed Model (PHMM)
#' 
#' Fits proportional hazards model incorporating random effects. The function
#' implements an EM agorithm using Markov Chain Monte Carlo at the E-step as
#' described in Vaida and Xu (2000).
#' 
#' \tabular{ll}{ Package: \tab phmm\cr Version: \tab 0.2\cr Date: \tab
#' 2008-01-15\cr Depends: \tab survival\cr Suggests: \tab lme4\cr License: \tab
#' GPL2\cr Packaged: \tab Fri Jul 11 10:33:57 2008; mdonohue\cr Built: \tab R
#' 2.8.0; universal-apple-darwin8.11.1; 2008-11-29 12:05:00; unix\cr }
#' 
#' Index: \preformatted{ AIC.phmm Akaike Information Criterion for PHMM cAIC
#' Conditional Akaike Information Criterion for PHMM e1582 Eastern Cooperative
#' Oncology Group (EST 1582) linear.predictors PHMM Design loglik.cond PHMM
#' conditional log-likelihood phmm Proportional Hazards Model with Mixed
#' Effects phmm-package Proportional Hazards Model with Mixed Effects
#' phmm.cond.loglik PHMM conditional log-likelihood phmm.design PHMM Design
#' pseudoPoisPHMM Pseudo poisson data for fitting PHMM via GLMM traceHat Trace
#' of the "hat" matrix from PHMM-MCEM fit }
#' 
#' @name phmm-package
#' @docType package
#' @author Ronghui Xu, Michael Donohue
#' 
#' Maintainer: Michael Donohue \email{mdonohue@@ucsd.edu}
#' @references Vaida, F. and Xu, R. "Proportional hazards model with random
#' effects", \emph{Statistics in Medicine,} 19:3309-3324, 2000.
#' 
#' Donohue, MC, Overholser, R, Xu, R, and Vaida, F (January 01, 2011).
#' Conditional Akaike information under generalized linear and proportional
#' hazards mixed models. \emph{Biometrika}, 98, 3, 685-700.
#' @keywords package
#' @seealso \code{\link{phmm}}
#' @useDynLib phmm
NULL



