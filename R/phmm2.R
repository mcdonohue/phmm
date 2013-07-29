phmm2 <- function (formula, data, subset, 
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
  #    print(temp)
  for (i in c('data', 'subset', 'na.action')){
    if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
  }
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
  #   print(flist)
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

  #print(m)
  W <- model.matrix(as.formula(paste("~", random[1])), m)[,,drop=F]
  cluster <- model.matrix(as.formula(paste("~", random[2])), m)[,-1,drop=F]
  if(ncol(cluster) != 1) stop("Unsupported cluster structure")
  nclust <- length(unique(cluster))
  n <- nrow(Z)
  nfixed <- ncol(Z)
  names.random <- colnames(W)

  nrandom <- ncol(W)
  if (nrandom ==0) stop("No random effects terms found")

  if(Sigma == "identity") Sigma0 <- diag(1, nrandom)
  invSigma <- solve(Sigma0)
  detSigma <- det(Sigma0)
  Sigma    <- solve(Sigma0)

  if(varcov == "diagonal"){ 
    varcov = 2
  }else{
    stop(paste("\nUnknown structure specified for var-covariance matrix:", 
    varcov))
  }

  if(verbose){
    cat("\nProportional Hazards Mixed-Effects Model fit with MCMC-EM\n")
    date0 <- date()
  }

  nobs <- n
  ncov <- nfixed
  nreff <- nrandom         
  x <- as.numeric(Y[,1])
  z <- as.matrix(Z, ncol=ncov)    #;print("z"); print(z)
  w <- as.matrix(W, ncol=nreff)           #;print(print"w"); print(w)
  delta <- as.integer(Y[,2])
  cluster <- as.integer(cluster)
  varcov <- as.integer(varcov)
  NINIT <- 10     

  emstep <- 0
  nobs <- n
  ncov <- nfixed
  nreff <- nrandom
  bhat <- matrix(0, nrow=nreff, ncol=nclust)
  sdbhat <- matrix(0, nrow=nreff, ncol=nclust)     
  steps <- matrix(0, nrow=MAXSTEP+1, ncol=(ncov+nreff))   	
  var <- matrix(0, nrow=(ncov+nreff+nobs), ncol=(ncov+nreff+nobs))
  verbose <- as.integer(verbose)
  bridgeC <- 1
  laplacetot <- 1
  llaplace <- 1
  limport <- 1
  lbridge <- 1
  lambda <- rep(0, n)
  Lambda <- rep(0, n)

  sum0 <- rep(0, (nobs+1))
  sum1 <- matrix(0, nrow=ncov, ncol=(nobs+1))
  sum2 <- array(list(rep(0, (nobs+1))), c(ncov, ncov))

  # x, delta, z, w are sorted by time   
  temp.x <- x
  temp.delta <- delta
  x <- x[order(temp.x, -temp.delta)]
  delta <- delta[order(temp.x, -temp.delta)]     
  cluster <- cluster[order(temp.x, -temp.delta)]
  z <- as.matrix(z[order(temp.x, -temp.delta),])
  w <- as.matrix(w[order(temp.x, -temp.delta),])         
  #print(cbind(x, delta, cluster, z, w))
  
  # xx, ddelta, cluster, rank, zz, ww are sorted by clusters     
  rank <- seq(1, nobs)
  xx <- x
  ddelta <- delta
  zz <- z
  ww <- w
  temp.cluster <- cluster
  xx <- xx[order(temp.cluster)]
  ddelta <- ddelta[order(temp.cluster)]
  cluster <- cluster[order(temp.cluster)]     
  rank <- rank[order(temp.cluster)]
  zz <- as.data.frame(zz[order(temp.cluster),])
  ww <- as.data.frame(ww[order(temp.cluster),])        

  #sort data first by time and -event, then by cluster     
  data <- as.data.frame(data[order(data$time, -data$event),])      

  clust_sizes <- as.data.frame(table(cluster))[,"Freq"]          
  # print("clust_sizes"); print(clust_sizes)
  clust_start <- cumsum(c(1, clust_sizes))               
  # print("clust_start"); print(clust_start)

  if (verbose) {cat("\nNo. of clusters: ", nclust)}         
  if (verbose) {cat(" \n Cluster sizes: ", clust_sizes, "\n")}    

  B <- matrix(0, nrow=nreff, ncol=nclust)     
  a <- matrix(0, nrow=nreff, ncol=nclust)     

  for (i in 1:nclust){
    for (j in 1:nreff){
      B[j, i] <- 0
      a[j, i] <- 0
      for(k in clust_start[i]:(clust_start[i] + clust_sizes[i] - 1)){
        a[j, i] <- a[j, i] + ddelta[k]*ww[k, j]
      }
    }
  }

  detconv <- rep(0, nclust)
  sumb <- matrix(0, nrow=nreff, ncol=nclust)     
  sumbb <- matrix(0, nrow=nreff, ncol=nclust)     
  eb <- matrix(0, nrow=nreff, ncol=nclust)
  vv <- matrix(0, nrow=nreff, ncol=nclust)     #in c code vv is noted as v
  
  sumv <- array(0, dim=list(nreff, nreff, nclust))
  condvar <- array(0, dim=list(nreff, nreff, nclust))
  invcondv <- array(0, dim=list(nreff, nreff, nclust))
  laplace <- rep(0, nclust) 
  
  # call EM()
  result.EM <- EM(ncov, nreff, Sigma, nobs, omega, betahat, sum0, sum1, z, delta, sum2, 
    lambda, Lambda, Lambexp, steps, NINIT, emstep, MAXSTEP, maxtime, nclust, 
    sumb, sumbb, B, clust_start, ww, Gbs, Gbsvar, a, rank, CONVERG, bhat, 
    sdbhat, sumv, invSigma, eb, b, condvar, detSigma, verbose, varcov, flist, data)

  betahat <- result.EM$betahat
  Lambexp <- result.EM$Lambexp
  lambda <- result.EM$lambda
  Sigma <- result.EM$Sigma
  invSigma <- solve(result.EM$Sigma)
  omega  <- result.EM$omega  
  B <- result.EM$B  

  print(betahat)

  # call Var()
  var.fit  <- .C("sVar",
    Gbs = as.integer(Gbs),
    nobs = as.integer(nobs),
    nreff = as.integer(nreff),
    ncov = as.integer(ncov),
    NINIT = as.integer(NINIT),
    betahat = as.double(betahat), #
    nclust = as.integer(nclust),
    sumb = as.vector(sumb),
    sumbb = as.vector(sumbb),
    B = as.vector(B),
    clust_start = as.integer(clust_start),
    wwv =  as.vector(as.matrix(ww)),
    Lambexp = as.vector(Lambexp),
    lambda = as.vector(lambda),
    av = as.vector(a),
    Sigma = as.vector(Sigma),
    invSigma = as.vector(invSigma),
    rank = as.double(rank), # 18
    zzv =  as.vector(as.matrix(zz)),
    delta =  as.integer(c(0, delta)),
    ddelta =  as.integer(c(0, ddelta)),
    omega =as.double(omega), # 22
    z = as.vector(as.matrix(z)),
    var = as.double(rep(0, (ncov+nreff+nobs+1)*(ncov+nreff+nobs+1))), 
    PACKAGE="phmm" )
  #print("Sigma")
  #print(Sigma)         
  #var <- var.fit$var
} # end phmm2()
  
  
  #######################################################################################
phmm.Estep.MCMC <- function(ncov, nreff, Sigma, nobs, NINIT, nclust, sumb, 
  sumbb, sumv, B, ww, Gbs, rank, CONVERG, clust_start, omega, Lambexp, vv, condvar, eb,
  detSigma, varcov, data, a)
{
  eb <- matrix(0, nrow=nreff, ncol=nclust)
  err <- integer(1)
  ninit <- 4
  dometrop <- 0
  npoint <- 100
  ncent <- 0
  neval <- integer(1)
  nsamp <- 1
  xinit <- rep(0, NINIT)
  xl <- -100
  xr <- 100
  xprev <- 0.0
  xsamp <- integer(1)
  temp <- integer(1)
  b <- rep(0, nreff)
  sum <- rep(0, nreff)
  alpha <- rep(0, nobs)
  oomega <- rep(0, nobs)
  myxinit <- matrix(0, nrow=nreff, ncol=ninit)
  vtemp <- matrix(0, nrow=nreff, ncol=nreff)
  vv <- matrix(0, nrow=nreff, ncol=nclust)

  for (i in 1:nclust){
    b <- B[ ,i]

    alpha[(clust_start[i]:(clust_start[i + 1] - 1))] <- 0
    #alpha[(clust_start[i]:(clust_start[i + 1] - 1))] <- b %*% (t(ww[(clust_start[i]:(clust_start[i + 1] - 1)),]))
    for (l in (clust_start[i]:(clust_start[i + 1] - 1))){
    alpha[l] <- 0
    for (d in 1:nreff){
    alpha[l] <- alpha[l] + b[d]*ww[l, d]  
    }
    }
    ww <- as.matrix(ww)
    invSigma <- Sigma

    fit  <- .C("phmm2",
      ii = as.integer(i),
      NINIT = as.integer(NINIT),
      Gbs = as.integer(Gbs),
      nreff = as.integer(nreff),
      nobs = as.integer(nobs),
      nclust = as.integer(nclust),
      alpha =  as.double(alpha),
      b = as.double(b),
      Lambexp = as.double(Lambexp),
      clust_start = as.double(clust_start),
      wwv =  as.double(ww),
      av = as.double(a),
      Sigma = as.double(Sigma),
      rank = as.double(rank),
      invSigma = as.double(invSigma),
      sum = as.double(sum),
      sumb = as.double(sumb),
      sumbb = as.double(sumbb),
      sumv = as.double(sumv),
      oomega =as.double(oomega), 
      PACKAGE="phmm")

    alpha <- fit$alpha
    b  <- fit$b
    Lambexp  <- fit$Lambexp
    sum  <- fit$sum
    sumb  <- matrix(fit$sumb, nrow=nreff, ncol=nclust)
    sumbb  <- matrix(fit$sumbb, nrow=nreff, ncol=nclust)
    sumv  <- array(fit$sumv, dim=list(nreff, nreff, nclust))    
    oomega  <- fit$oomega    

    B[, i] <- b # B keeps current set of r.effs

    eb[, i] <- sumbb[, i]/Gbs
    vv[, i] <- sumb[, i]/Gbs - (sumbb[, i]/Gbs)*(sumbb[, i]/Gbs)
    condvar[,,i] <- sumv[,,i]/Gbs - (sumbb[ , i]/Gbs)%*%t(sumbb[ , i]/Gbs)
  } #end i

  
  omega[rank[1:nobs]] <- oomega[1:nobs]/Gbs
  temp <- array(0, dim=list(nreff, nreff))
  for (i in 1:nclust){temp <- temp + sumv[,,i]}
  Sigma <- temp/(nclust*Gbs)
  for (d in 1:nreff) {
    for (dd in 1:nreff){
      if (d!=dd) Sigma[d,dd]  <- 0
    }
  }
  vtemp <- Sigma
  detSigma <- det(vtemp) 
  list(omega=omega, Sigma=Sigma, B=B)
}

  
####################################################################################
#EM()  
EM <- function(ncov, nreff, Sigma, nobs, omega, betahat, sum0, sum1, z, delta, sum2, 
  lambda, Lambda, Lambexp, steps, NINIT, emstep, MAXSTEP, maxtime, nclust, 
  sumb, sumbb, B, clust_start, ww, Gbs, Gbsvar, a, rank, CONVERG, bhat, 
  sdbhat, sumv, invSigma, eb, b, condvar, detSigma, verbose, varcov, flist, data)
{
  time_start <- Sys.time()  
  beta <- rep(0, ncov)
  temp_G <- matrix(0, nrow=nreff, ncol=nreff)     
  omega <- rep(1, nobs)
  betahat <- rep(0, ncov)     
  fmla <- flist$fixed

  start.coxph <- coxph(fmla, data, model=TRUE)

  betahat <- start.coxph$coefficient
  Lambda <- basehaz(start.coxph, centered=FALSE)$hazard
      
  Lambda <- rep(Lambda, table(data$time))  # handle ties 
  lambda <- c(Lambda[1], diff(Lambda))
  aa <- rep(0, nobs) 
  aa <- z%*%betahat  
  
  Lambexp <- Lambda*exp(aa) 

  steps[emstep+1, c(1:ncov)] <- betahat
  steps[emstep+1, c((ncov+1):(ncov+nreff))] <- diag(Sigma)

  if (verbose) cat("\nstep", paste("beta", 1:ncov, sep=""), paste("var", 1:nreff, sep=""), ": \n", sep="   ")

  if (verbose) cat(" ", 0, sprintf("%.4f", betahat), sprintf("%.4f", diag(Sigma)), "\n", sep="  ")

  # Test for EM convergence
  repeat {
    emstep <- emstep + 1
    if(emstep > MAXSTEP) {
      for (i in 1:nclust){
        for (d in 1:nreff){
          bhat[d, i]   <- sumbb[d, i]/Gbs
          sdbhat[d, i] <- sqrt(sumb[d, i]/Gbs - (sumbb[d, i]/Gbs)*(sumbb[d, i]/Gbs))
        }
      } 
      break
    }
  
    #Gibbs E-step
    #print(a)
    #print("phmm.Estep.MCMC calleds")
    #t.start <- Sys.time()

    #print("####################################")
    #print("ncov");   print(ncov)
    #print("nreff");       print(nreff)
    #print("Sigma");   print(Sigma)
    #print("nobs");     print(nobs)
    #print("NINIT");   print(NINIT)
    #print("nclust");    print(nclust)
    #print("Lambexp");      print(Lambexp)
    #print("clust_start");  print(clust_start)
    #print("ww");           print(ww)
    #print("a");    print(a)
    #print("Sigma");        print(Sigma)
    #print("rank");         print(rank)
    #print("invSigma");     print(invSigma)
    #print("sum");          print(sum)
    #print("sumb");         print(sumb)
    #print("sumbb");        print(sumbb)
    #print("sumv");         print(sumv)
    #print("oomega");       print(oomega)    

    erlt <- phmm.Estep.MCMC(ncov, nreff, Sigma, nobs, NINIT, nclust, sumb, 
      sumbb, sumv, B, ww, Gbs, rank, CONVERG, clust_start,  omega, Lambexp, vv, condvar, 
      detSigma, varcov, data, a=a)
    #t.stop <- Sys.time()
    #t.time <- t.stop-t.start
    #print(t.time)
    #print("Sigma after e step") ; print(Sigma)
    #M-step for beta
    oomega  <- erlt$oomega
    omega  <- erlt$omega
    # print("omega")
    # print(omega)
    fmla.fit.coxph <- update(fmla, ~. + offset(log.omega))
    log.omega <- log(erlt$omega)      
    # print("log.omega") 
    # print(log.omega)
    Sigma <- erlt$Sigma        
    #print(Sigma)
    B <- erlt$B       
    dataset <- data       
    dataset <- cbind(dataset, log.omega)
    fit.coxph <- coxph(fmla.fit.coxph, dataset, model=TRUE, method=c("breslow"))

    betahat <- fit.coxph$coefficient
    Lambda <- basehaz(fit.coxph, centered=FALSE)$hazard           
    Lambda <- rep(Lambda, table(data$time))
    # print("Lambda")
    # print(Lambda)  
    lambda <- c(Lambda[1], diff(Lambda))  
    # print("lambda")
    # print(lambda)  
    aa <- rep(0, nobs) 
    aa <- z%*%betahat
    Lambexp <- Lambda*exp(aa)

    steps[emstep+1, c(1:ncov)] <- betahat
    steps[emstep+1, c((ncov+1):(ncov+nreff))] <- diag(erlt$Sigma)  # ;print(steps)

    if (verbose) cat(format(emstep, trim=TRUE, justify="right", width=4), sprintf("%.4f", betahat), sprintf("%.4f", diag(erlt$Sigma)), "\n", sep="  ")  

    time_stop <- Sys.time()       
    dif <- difftime(time_stop, time_start)

    if (dif > maxtime){
      emstep <- CONVERG
      CONVERG <- CONVERG + 1
    }

    if (emstep > CONVERG) {Gbs <- Gbsvar}     

  }# end repeat
  return(list(betahat=betahat, Lambexp=Lambexp, lambda=lambda, omega=omega, Sigma=Sigma, B=B, oomega=oomega))
}# end EM
  



