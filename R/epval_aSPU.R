epval_aSPU <- function(sam1, sam2, pow = c(1:6, Inf), perm.iter = 1000){ 
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	pval <- aSPUperm(sam1, sam2, pow = pow, n.perm = perm.iter)
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pow <- pow
	out$pval <- pval
	return(out)
}