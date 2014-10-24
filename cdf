Lcdf <- function(x, baseline=NULL, q=c(0,0.5,0.95)) {	
		
	if(is.null(baseline)) {
		baseline <- min(x)
		n <- length(x) -1
		}
	else n <- length(x)
	
	X <- x - baseline
	
	L <- max(X)/(1-q)^(1/n) + baseline
	
	names(L) <- paste(q*100, "%", sep='')
	return(L)
}
