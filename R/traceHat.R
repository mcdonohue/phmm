cAIC <- function(object, ..., k = 2) UseMethod("cAIC")

cAIC.phmm <- function(object, method = "direct", ..., k = 2){
  as.numeric(-2*object$loglik["Conditional"])+k*traceHat(object, method = method)
}

cAIC.coxph <- function(object, ..., k = 2){
	-2*object$loglik[2] + k*length(object$coef)
}

AIC.coxph <- cAIC.coxph

traceHat <- function(x, method = "direct"){ 
  if(! method %in% c("direct", "pseudoPois", "HaLee", "approx")) 
    stop("Undefined traceHat method.")

    time <- x$Y[, 1]; delta <- x$Y[, 2]; z <- x$Z
    w <- x$W; b <- x$bhat.long; Sigma <- Matrix(x$Sigma)
    eventtimes <- unique(sort(time[delta == 1]))
    neventtimes <- length(eventtimes)
    cluster <- as.numeric(as.character(x$cluster))
    nclust <- length(unique(cluster))

  if(method%in%c("direct", "HaLee", "approx")){
    # use Ha, Lee, MacKenzie 2007 framework
    # to compute traceHat directly
    # z -- covs for fixed effects
    # w -- covs for random effects
    # time -- min(T, C)
    # Sigma -- varcov Matrix for random effects
    # fitted -- linear part

    fitted <- x$linear.predictors  
    xx <- cBind(ID = 1:length(time), time = time, delta = delta)
    Lambda = cumsum(x$lambda)[!duplicated(sort(time), fromLast = TRUE)]
    lambda = Lambda-c(0, Lambda[1:(length(Lambda)-1)])
    Lambda <- cBind(Lambda = Lambda, 
            lambda = lambda, 
                    time = unique(sort(time)))
    xx <- merge(xx, Lambda, by = "time", all.x = TRUE)
    Lambda <- xx[order(xx$ID), "Lambda"]
    lambda <- xx[order(xx$ID), "lambda"]
  
    xxx <- xx[xx$delta == 1, c("time", "lambda")]
    ulambda <- xxx[!duplicated(xxx$time), "lambda"]
  
    time <- time[order(cluster)]
    delta <- delta[order(cluster)]
    if(ncol(z) == 1){ z <- Matrix(z[order(cluster)])
      }else{ z <- Matrix(z[order(cluster), ]) }
    if(ncol(w) == 1){ 
      w <- Matrix(w[order(cluster)])
      b <- Matrix(b[order(cluster)])
      bhat <- b[!duplicated(sort(cluster))]
      }else{ 
        w <- Matrix(w[order(cluster), ]) 
        b <- Matrix(b[order(cluster), ])
        bhat <- b[!duplicated(sort(cluster)), ]
        }
    
    cluster <- sort(cluster)
    nclust <- length(unique(cluster))
    cluster <- rep(1:nclust, table(cluster))
  
    X <- z
    Z <- bdiag(lapply(unique(cluster), function(x){
         block <- w[cluster == x, ]
         if(ncol(w) > 1) return(rBind(block)) else return(cBind(block))
       }))
    D <- bdiag(rep(list(Sigma), nclust))
  
    W3 <- diag(as.vector(exp(fitted)))
    B <- diag(Lambda)
    W1 <- W3 %*% B
    
    M <- outer(time, eventtimes, FUN = function(x, y) ifelse(x >= y, 1, 0))
  
    dk <- rep(0, neventtimes)
    for(k in 1:neventtimes)
      dk[k] <- sum(delta[time == eventtimes[k]])
    C <- diag(dk/(ulambda^2))
    W2 <- (W3 %*% M) %*% solve(C) %*% t(W3 %*% M)
    W <- W1-W2
    
    if(nrow(Sigma) == 1 & Sigma[1, 1] == 0){ 
      return(sum(diag(solve(crossprod(X, W) %*% X %*% crossprod(X, W) %*% X))))
    }else if(method == "direct"){
      U <- cBind(as.matrix(X), as.matrix(Z))
      A <- Matrix(bdiag(matrix(0, x$nfixed, x$nfixed), solve(D)))
      UW <- crossprod(U, W)
      return(sum(diag(UW %*% U %*% 
                solve(UW %*% U + A))))
    }else if(method %in% c("HaLee", "approx")){
      U <- solve(D)
      XW <- crossprod(X, W); ZW <- crossprod(Z, W)
      J <- if(method == "HaLee") 0 else U %*% as.vector(t(bhat)) %*% t(as.vector(t(bhat))) %*% U
      return(sum(diag(solve(
      rBind(cBind(XW %*% X, XW %*% Z), 
            cBind(ZW %*% X, ZW %*% Z + U))) %*% 
      rBind(cBind(XW %*% X, XW %*% Z), 
            cBind(ZW %*% X, ZW %*% Z + J)))))
    }
  }else if(method == "pseudoPois"){
    xx <- pseudoPoisPHMM(x)
    
    xx <- cBind( ID = 1:nrow(xx), xx)
    xx <- xx[order(xx[, "cluster"], xx[, "ID"]), ]
    z <- xx[, c(paste("t", 1:neventtimes, sep = ''), 
                paste("z", 1:x$nfixed, sep = ''))]
    w <- Matrix(as.matrix(xx[, paste("w", 1:x$nrandom, sep = '')]))
    cluster <- xx[, "cluster"]
    fitted <- xx[, "linear.predictors"] + 
      xx[, paste("t", 1:neventtimes, sep = '')] %*% log(x$lambda[x$lambda != 0]) + 
      log(xx[, "N"])
    nclust <- length(unique(cluster))
    cluster <- rep(1:nclust, table(cluster))
    Sigma = Matrix(x$Sigma)

    X <- as.matrix(z)
    Z <- bdiag(lapply(unique(cluster), function(x) w[cluster == x, ]))
    D <- bdiag(rep(list(Sigma), nclust))
    W <- Diagonal(x = as.numeric(exp(fitted)))
    U <- solve(D)
    ZWZ <- crossprod(Z, W) %*% Z
    ZWX <- crossprod(Z, W) %*% X
    XWX <- crossprod(X, W) %*% X
    XWZ <- crossprod(X, W) %*% Z

    return(nclust * x$nrandom + x$nfixed - sum(diag(solve(ZWZ + U - ZWX %*% solve(XWX) %*% XWZ) %*% U)))
    }else return(NULL) # end pseudoPois
}
