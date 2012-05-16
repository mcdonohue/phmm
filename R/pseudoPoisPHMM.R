# creates pseudo data for poisson PHMM model
# from object of class "phmm"
pseudoPoisPHMM <- function (x) UseMethod("pseudoPoisPHMM")
pseudoPoisPHMM.phmm <- function(x){	
	dd <- cBind(x$cluster, x$Z, x$W)
	group <- apply(dd,1,paste,collapse="XX")
	groups <- unique(group)
	dd <- cBind(x$Y[, 1], x$Y[, 2], x$linear.predictors, dd)
	colnames(dd) <- c("time", "delta", "linear.predictors",
		"cluster",
		paste("z", 1:x$nfixed, sep=''),
		paste("w", 1:x$nrandom, sep=''))
	ddext <- lapply(sort(unique(dd[dd[,"delta"]==1,"time"])), function(t){
		tdd <- lapply(1:length(groups), function(i){
			unlist(c(t,sum(dd[group==groups[i],"time"]>=t),
			sum(dd[group==groups[i]&dd[,"delta"]==1,"time"]==t),
			dd[group==groups[i]&!duplicated(group==groups[i]),
				c("cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")]))
			})
		})
	cnames <- c('time','N','m',"cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")
	ddext <- Matrix(unlist(ddext), byrow=TRUE, ncol=length(cnames), sparse=TRUE)
	colnames(ddext) <- cnames
	times <- sort(unique(dd[dd[,"delta"]==1,"time"]))
	timematrix <- Matrix(0,nrow(ddext),ncol=length(times))
	colnames(timematrix) <- paste("t",1:length(times),sep='')
	for(i in 1:length(times) ){		
		timematrix[ddext[,"time"]==times[i],paste("t",i,sep='')] <- 1
	}
	ddext <- cBind(ddext, timematrix)
	ddext <- ddext[ddext[,'N']!=0,]
	ddext <- ddext[order(ddext[,'cluster']),]
	return(ddext)
}

pseudoPoisPHMM.coxph <- function(x){	
	group <- apply(x$x,1,paste,collapse="XX")
	groups <- unique(group)
	dd <- cBind(x$y[, 1], x$y[, 2], x$x%*%x$coef, x$x)
	colnames(dd) <- c("time", "delta", "linear.predictors",
		paste("z", 1:ncol(x$x), sep=''))
	ddext <- lapply(sort(unique(dd[dd[,"delta"]==1,"time"])), function(t){
		tdd <- lapply(1:length(groups), function(i){
			unlist(c(t,sum(dd[group==groups[i],"time"]>=t),
			sum(dd[group==groups[i]&dd[,"delta"]==1,"time"]==t),
			dd[group==groups[i]&!duplicated(group==groups[i]),
				c(paste("z", 1:ncol(x$x), sep=''),
				  "linear.predictors")]))
			})
		})
	cnames <- c('time','N','m',
				paste("z", 1:ncol(x$x), sep=''),
				"linear.predictors")
	ddext <- Matrix(unlist(ddext), byrow=TRUE, ncol=length(cnames), sparse=TRUE)
	colnames(ddext) <- cnames
	times <- sort(unique(dd[dd[,"delta"]==1,"time"]))
	timematrix <- Matrix(0,nrow(ddext),ncol=length(times))
	colnames(timematrix) <- paste("t",1:length(times),sep='')
	for(i in 1:length(times) ){		
		timematrix[ddext[,"time"]==times[i],paste("t",i,sep='')] <- 1
	}
	ddext <- cBind(ddext, timematrix)
	ddext <- ddext[ddext[,'N']!=0,]
	# bh <- basehaz(x, centered = FALSE)
	# lambda <- bh$hazard - c(0,bh$hazard[1:(length(bh$hazard)-1)])
	# alpha <- log(ddext[,paste("t",1:length(times),sep='')]%*%lambda)
	# ddext <- cBind(ddext, alpha)
	return(ddext)
}