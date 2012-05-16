phmm <- function (formula, data, subset, 
	na.action = na.fail, Sigma = "identity", varcov = "diagonal", 
	NINIT = 10, VARSTART = 1, MAXSTEP = 100, CONVERG = 90, Gbs = 100, 
	Gbsvar = 1000, verbose = FALSE, maxtime = 120, random)
{
    Call <- match.call()

    if (!missing(random)) {
        warning("The 'random' argument of phmm is deprecated")
        if (class(random) != 'formula' || length(random) !=2) 
            stop("Invalid random formula")
        j <- length(formula)   #will be 2 or 3, depending on if there is a y

        # Add parens to the random formula and paste it on
        formula[[j]] <- call('+', formula[[j]], call('(', random[[2]]))  
        }

    temp <- call('model.frame', formula= subbar(formula))
    for (i in c('data', 'subset', 'na.action'))
        if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
        Y <- model.extract(m, "response")
        n <- nrow(Y)
        if (!inherits(Y, "Surv")) stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type!='right')
            stop(paste("phmm doesn't support '", type,
                              "' survival data", sep=''))

        # Check for penalized terms; the most likely is pspline
        pterms <- sapply(m, inherits, 'coxph.penalty')
        if (any(pterms)) {
            stop("Penalized terms are not supported by phmm")
            }

        flist <- formula1(formula)
        if (hasAbar(flist$fixed))
            stop("Invalid formula: a '|' outside of a valid random effects term")

        special <- c("strata", "cluster")
        Terms <- terms(flist$fixed, special)
        attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
        strats <- attr(Terms, "specials")$strata
        cluster<- attr(Terms, "specials")$cluster
        if (length(cluster)) {
            stop ("A cluster() statement is invalid in phmm")
            }
        if (length(strats)) {
            stop ("A strata() statement is not supported in phmm")
            }
        else Z <- model.matrix(Terms, m)[,-1,drop=F]
    
    
    random <- unlist(strsplit(as.character(flist$random[[1]])[2], "|", fixed=TRUE))
    W <- model.matrix(as.formula(paste("~", random[1])), m)[,,drop=F]
    cluster <- model.matrix(as.formula(paste("~", random[2])), m)[,-1,drop=F]
    if(ncol(cluster) != 1) stop("Unsupported cluster structure")
    nclust <- length(unique(cluster))
    n <- nrow(Z)
    nfixed <- ncol(Z)
	names.random <- colnames(W)

    nrandom <- ncol(W)
    if (nrandom ==0) stop("No random effects terms found")

	if(Sigma == "identity") Sigma0 = diag(1, nrandom)
	invSigma = rbind(0, cbind(0, solve(Sigma0)))
	detSigma = det(Sigma0)
	Sigma    = rbind(0, cbind(0, solve(Sigma0)))

	if(varcov == "diagonal"){ varcov = 2
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:", 
		varcov))}

    if(verbose){
		cat("\nProportional Hazards Mixed-Effects Model fit with MCMC-EM\n")
		date0 <- date()
	}

	fit <- .C("phmm", 
		X = as.numeric(c(0,Y[,1])), 
		Z = rbind(0, as.matrix(Z)), 
		W = rbind(0, as.matrix(W)),
		delta = as.integer(c(0,Y[,2])), 
		cluster = as.integer(c(0,cluster)),
		varcov = as.integer(varcov),
		Sigma = Sigma,
		invSigma = invSigma,
		detSigma = as.double(detSigma),
		NINIT=as.integer(10),
		MAXSTEP=as.integer(MAXSTEP),
		maxtime=as.double(maxtime),
		CONVERG=as.integer(CONVERG),
		emstep = as.integer(0),
		Gbs = as.integer(Gbs),
		Gbsvar = as.integer(Gbsvar),
		n = as.integer(n),
		nfixed = as.integer(nfixed),
		nrandom = as.integer(nrandom),
		bhat = double(nclust*nrandom),
		sdbhat = double(nclust*nrandom),
		steps = double((MAXSTEP+1)*(nfixed+nrandom)),
		var = double((nfixed+nrandom+n+1)^2),
		verbose = as.integer(verbose),
		bridgeC = double(1),
		laplacetot = double(1),
		llaplace = double(1), 
		limport = double(1), 
		lbridge = double(1),
		lambda = double(n+1),
		Lambda = double(n+1),
		PACKAGE="phmm" )

	if(varcov == 2) fit$varcov = "diagonal"
	
	class(fit) <- "phmm"

	fit$Call <- Call
	fit$formula <- formula
	fit$random <- random
	fit$verbose <- as.logical(fit$verbose)
	fit$steps <- matrix(fit$steps, nrow = MAXSTEP+1, byrow = TRUE)
	colnames(fit$steps) <- c(colnames(Z), paste("var(", names.random, ")", sep = ""))
	rownames(fit$steps) <- 0:MAXSTEP
	
	fit$var <- matrix(fit$var, nrow = nfixed+nrandom+n+1, byrow = FALSE)
	fit$var <- fit$var[1 + 1:(nfixed+nrandom), 1 + 1:(nfixed+nrandom)]
	rownames(fit$var) <- colnames(fit$var) <- c(colnames(Z), names.random)
	fit$varFix <- fit$var[1:nfixed, 1:nfixed]
	
	fit$Sigma0 <- Sigma0
	fit$Sigma <- matrix(fit$Sigma, nrow = nrandom+1, byrow = FALSE)
	fit$Sigma <- fit$Sigma[1 + 1:nrandom, 1 + 1:nrandom]
	if(nrandom == 1){ names(fit$Sigma) = paste("var(", names.random, ")", sep = "")
	}else{
		colnames(fit$Sigma) = names.random
		rownames(fit$Sigma) = names.random
	}
	fit$invSigma <- matrix(fit$invSigma, nrow = nrandom+1, byrow = FALSE)
	fit$invSigma <- fit$invSigma[1 + 1:nrandom, 1 + 1:nrandom]

	fit$coefficients = fit$steps[MAXSTEP+1, 1:fit$nfixed]
	names(fit$coefficients) <- colnames(Z)
	fit$bhat <- matrix(fit$bhat, nrow = nclust, byrow = TRUE)
	fit$sdbhat <- matrix(fit$sdbhat, nrow = nclust, byrow = TRUE)
	
	rownames(fit$bhat) <- rownames(fit$sdbhat) <- unique(cluster) 
	colnames(fit$bhat) <- colnames(fit$sdbhat) <- colnames(W)	
	
	fit$X = Y[, 1]
	fit$Z = as.matrix(Z)
	fit$W = as.matrix(W)
	fit$delta = Y[, 2] 
	fit$cluster = cluster
	
	fit$Lambda = fit$Lambda[-1]
	fit$lambda = fit$lambda[-1]
	
	fit$bhat.long <- as.matrix(merge(
	   cbind(cluster = as.numeric(rownames(fit$bhat)), as.matrix(fit$bhat)), 
	   cbind(cluster = as.numeric(as.character(fit$cluster)))))[, -1]
	fit$Y <- Y
	
	if(fit$nrandom > 2) fit$lbridge <- NULL
	fit$loglik <- c(Conditional = loglik.cond(fit), 
		Laplace = fit$llaplace, RIS = fit$limport, BS = fit$lbridge)
		
	fit$linear.predictors <- linear.predictors(fit)
	fit$flist <- flist
	fit$Terms <- Terms

	fit <- fit[c('Call', 'Y', 'Z', 'W', 'cluster', 'varcov', 'Sigma', 
	   'invSigma', 'detSigma', 'NINIT', 'MAXSTEP', 'CONVERG', 'emstep', 
	   'Gbs', 'Gbsvar', 'n', 'nfixed', 'nrandom', 'bhat', 'bhat.long', 
	   'sdbhat', 'steps', 'var', 'verbose', 'loglik', 'lambda', 'Lambda', 
	   'varFix', 'coefficients', 'linear.predictors',
	   'formula', 'flist', 'Terms')]
	
	class(fit) <- "phmm"
    return(fit)
}

loglik.cond <- function (x) UseMethod("loglik.cond")
loglik.cond.phmm <- function(x){
	#Function to compute conditional log-likelihood
	phmm.cond.loglik(time = x$Y[, 1], delta = x$Y[, 2], z = x$Z, beta = x$coef, w = x$W, b = as.matrix(x$bhat.long))
}

phmm.cond.loglik <- function(time, delta, z, beta, w, b){
	#Function to compute conditional log-likelihood
    z <- as.matrix(z)
    wb <- matrix(0, nrow = nrow(w), ncol = ncol(w))
    for(i in 1:ncol(w)){
	  if(length(w[, i]) != length(b[, i])) stop("length(w[, i]) != length(b[, i])")
	  wb[, i] <- w[, i]*b[, i]
	  }
    if(!is.null(dim(wb))) wb <- apply(wb, 1, sum)
    numerator <- exp(z%*%beta+wb)
    denominator <- unlist(lapply(time, 
      FUN = function(x){
          sum(exp((z%*%beta+wb)[time >= x]))
          }))
    sum(ifelse(delta, 1, 0)*log(numerator/denominator))
}

AIC.phmm <- function(object, ..., k = 2){
	if(object$varcov == "diagonal"){ 
		return(-2*object$loglik+k*(object$nrandom+object$nfixed))
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:", 
		object$varcov))}
}
	
print.phmm <-
 function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nProportional Hazards Mixed-Effects Model fit by MCMC-EM\n")
  cat("  Model:", deparse(x$formula), "\n")
  cat("  Data:", deparse( x$Call$data ), "\n")
  cat("  Log-likelihood:\n")
	print(x$loglik[!is.null(x$loglik)], digits = digits, ...)
  cat("\nFixed effects:", deparse(x$flist$fixed), "\n")
  print(coef(x), digits = digits, ...)
  cat("\n")
  cat("Random effects:", deparse(x$flist$random[[1]]), "\n")
  cat("Estimated variance-covariance matrix:\n")
  print(x$Sigma, digits = digits, ...)
#  cat("Variance-Covariance:\n")
#	print(x$var, ...)
  cat("\nNumber of Observations:", x$n)
  cat("\nNumber of Groups: ", nrow(x$bhat))
  cat("\n\n")
}

print.summary.phmm <- print.phmm

summary.phmm <-
 function(object, ...)
{
	object$coefficients = cbind(Estimate = object$coef, 
	   Std.Error = sqrt(ifelse(diag(object$var)<0, NA, diag(object$var)))[1:object$nfixed])
	class(object)<-"summary.phmm"
	return(object)
}

plot.phmm <-
 function(x, ...)
{
	x = as.data.frame(x$steps)
	colnames(x) = make.names(colnames(x))
	fm = paste(paste(colnames(x), collapse = ' + '), "EM.Step", sep = " ~ ")
	x$EM.Step = as.numeric(rownames(x))
	xyplot(formula(fm), data = x, type = "l", allow.multiple = TRUE, outer = TRUE, scales = "free", ...)
}