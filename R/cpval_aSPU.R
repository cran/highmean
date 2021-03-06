cpval_aSPU <- function(sam1, sam2, pow = c(1:6, Inf), n.iter = 1000, seeds){
	if(missing(seeds)){
		seeds <- NULL
	}else{
		if(length(seeds) != n.iter){
			seeds <- NULL
			cat("The length of seeds does not match the specified n.iter.\n")
			cat("Seeds for each permutation/resampling iteration are assigned randomly.\n")
		}
	}

	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1*n2)/(n1 + n2)
	p <- dim(sam1)[2]

	if(!is.element(Inf, pow)) pow <- c(pow, Inf)
	pow <- sort(pow)

	diff <- colMeans(sam1) - colMeans(sam2)
	f.pow <- pow[pow != Inf]
	sam <- rbind(sam1, sam2)
	Ts.perm <- matrix(NA, length(pow), n.iter)
	for(b in 1:n.iter){
		if(!is.null(seeds)) set.seed(seeds[b])
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		diff.perm <- colMeans(sam1.perm) - colMeans(sam2.perm)
		for(j in 1:length(f.pow)){
			Ts.perm[j, b] <- sum(diff.perm^pow[j])
		}
	}

	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	pval <- numeric(length(pow) + 1)
	L <- numeric(length(pow))
	L.e <- numeric(length(pow) - 1)
	L.var <- numeric(length(pow) - 1)
	stan.L <- numeric(length(pow) - 1)
	for(i in 1:length(pow)){
		if(pow[i] != Inf){
			L[i] <- sum(diff^(pow[i]))
			L.e[i] <- mean(Ts.perm[i,])
			L.var[i] <- var(Ts.perm[i,])
			stan.L[i] <- (L[i] - L.e[i])/sqrt(L.var[i])
			if(pow[i] %% 2 == 1) pval[i] <- 2*(1 - pnorm(abs(stan.L[i])))
			if(pow[i] %% 2 == 0) pval[i] <- 1 - pnorm(stan.L[i])
		}else{
			diag.sam.cov <- diag(sam.cov)
			diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
			L[i] <- max(diff^2/diag.sam.cov)
			stan.L.inf <- tau*L[i] - (2*log(p) - log(log(p)))
			pval[i] <- pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(3.14159))
		}
	}

	odd.ga <- f.pow[f.pow %% 2 == 1]
	odd.ga.id <- which(f.pow %% 2 == 1)
	even.ga <- f.pow[f.pow %% 2 == 0]
	even.ga.id <- which(f.pow %% 2 == 0)
	n.odd.ga <- length(odd.ga)
	n.even.ga <- length(even.ga)
	R_O <- matrix(NA, n.odd.ga, n.odd.ga)
	R_E <- matrix(NA, n.even.ga, n.even.ga)
	diag(R_O) <- 1
	diag(R_E) <- 1
	for(s in odd.ga){
		for(t in odd.ga){
			if(s != t){
				io <- which(odd.ga == s)
				jo <- which(odd.ga == t)
				i <- which(pow == s)
				j <- which(pow == t)
				L.cov <- cov(Ts.perm[i,], Ts.perm[j,])
				R_O[io, jo] <- L.cov/sqrt(L.var[i]*L.var[j])
			}
		}
	}
	for(s in even.ga){
		for(t in even.ga){
			if(s != t){
				ie <- which(even.ga == s)
				je <- which(even.ga == t)
				i <- which(pow == s)
				j <- which(pow == t)
				L.cov <- cov(Ts.perm[i,], Ts.perm[j,])
				R_E[ie, je] <- L.cov/sqrt(L.var[i]*L.var[j])
			}
		}
	}
	TO <- max(abs(stan.L[odd.ga.id]))
	TE <- max(stan.L[even.ga.id])
	pval_O <- 1 - pmvnorm(lower = -rep(TO, n.odd.ga), upper = rep(TO, n.odd.ga),
		mean = rep(0, n.odd.ga), sigma = R_O)
	pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(TE, n.even.ga),
		mean = rep(0, n.even.ga), sigma = R_E)
	pval.min <- min(c(pval_O, pval_E, pval.inf))
	pval[length(pow) + 1] <- 1 - (1 - pval.min)^3
	names(pval) <- c(paste("SPU_", pow, sep = ""), "aSPU")

	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pow <- pow
	names(L) <- pow
	out$spu.stat <- L
	names(L.e) <- names(L.var) <- pow[pow < Inf]
	out$spu.e <- L.e
	out$spu.var <- L.var
	colnames(R_O) <- rownames(R_O) <- odd.ga
	out$spu.corr.odd <- R_O
	colnames(R_E) <- rownames(R_E) <- even.ga
	out$spu.corr.even <- R_E
	out$cov.assumption <- "the two groups have same covariance"
	out$method <- "combination of permutation method and asymptotic distributions"
	out$pval <- pval
	return(out)
}