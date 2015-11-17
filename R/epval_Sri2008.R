epval_Sri2008 <- function(sam1, sam2, perm.iter = 1000){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	tau <- (n1 + n2)/(n1*n2)
	n <- n1 + n2 - 2
	p <- dim(sam1)[2]
	Sn <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	diag.s <- diag(Sn)
	diag.s[diag.s == 0] <- 10^(-10)
	sam <- rbind(sam1, sam2)
	diff <- colMeans(sam1) - colMeans(sam2)
	XDX <- sum(diff^2/diag.s)
	sri2008.stat <- as.numeric(XDX)

	sri2008.stat.perm <- numeric(perm.iter)
	for(i in 1:perm.iter){
		perm <- sample(1:(n1 + n2))
		sam.perm <- sam[perm,]
		sam1.perm <- sam.perm[1:n1,]
		sam2.perm <- sam.perm[(n1 + 1):(n1 + n2),]
		Sn.perm <- ((n1 - 1)*cov(sam1.perm) + (n2 - 1)*cov(sam2.perm))/(n1 + n2 - 2)
		diag.s.perm <- diag(Sn.perm)
		diag.s.perm[diag.s.perm == 0] <- 10^(-10)
		diff.perm <- colMeans(sam1.perm) - colMeans(sam2.perm)
		XDX.perm <- sum(diff.perm^2/diag.s.perm)
		sri2008.stat.perm[i] <- as.numeric(XDX.perm)
	}
	sri2008.pval <- sum(sri2008.stat.perm >= sri2008.stat)/perm.iter
	names(sri2008.pval) <- "Sri2008"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pval <- sri2008.pval
	return(out)
}