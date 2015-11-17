apval_Chen2014 <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)
	test.stat <- max_threshold(sam1, sam2, sam.cov)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - exp(-exp(-test.stat))
	names(pval) <- "Chen2014"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$pval <- pval
	return(out)
}