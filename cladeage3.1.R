##################################################################
### FUNCTIONS FOR ESTIMATING CLADE AGE USING THE FOSSIL RECORD ###
##################################################################

### Author ###
# Santiago Claramunt
# Department of Natural History, Royal Ontario Museum, Toronto, Ontario M5S 2C6, Canada.
# E-mail: sclaramunt@rom.on.ca

#Language: R (https://www.r-project.org)


### Details ###

# These functions model the ucertainty about the age of a clade based on its fossil record. The primary strategy is to generate a random sample of possible clade ages using quantile functions. These numbers can then be used to generate a histogram and fit probability functions for using as priors in Bayesian time tree estiamtion (Claramunt & Cracraft 2015).


### Arguments ###

# ages		vector of fossil ages (age.lik) or a matrix (PD.age.unc) of fossil ages in which the first column is the maximum age (or lower stratigraphic bound) and the second column is the minimum fossil age (upper stratigraphic bound) for each fossil. 

# baseline 	baseline like the present (=0), or some other age of reference. If a baseline is not provided, the minimum value of age is used as baseline n-1 is used instead of n (Solow 2003)

# rescale	If TRUE (default) rescale the likelihood values to maje the curves comprarable across pseudoreplicates with slightly different ages.

# reps		number of pseudoreplicates used to account for fossil age uncertainty.

# line.col	color of the line plot (PD.age.plot)

# ...		further arguments passed to plot()




### FUNCTIONS ###

### Quantile function that computes the age correspoding to a particular probability for the upper bound of a distribution of ages (Strauss & Sadler 1989, Gingerich & Uhen 1998, Solow 2003). The method of Strauss & Sadler (1989) asumes that the distribution of fossil ages is uniform and their formula depends on the fossil ages range and the number of fossil ages. The method of Solow (2003) is a general method for non-uniform distributions and depends on the temporal gap between the oldest and the second oldest fossil ages. Both methods assume that fossil ages are independent samples from the same distribution (only relevant for the two oldest ages for Solow's method), therefore, fossils should be as independent as possible (ideally from different geological formations and different regions).

# In the particular case where there are only two fossil ages and no baseline age specified, Strauss & Sadler's and Solow's methods converge to the same result; the quantile functions are simply Xn/(1-P), and the likelihood function is 1/X.


# Argumetns #

# p		The desired probability level

# ages	Vector of fossil ages.

# baseline 	Youngest bound of the distribution that could be the present (=0) or some other geological age of reference. If a baseline is not provided, the minimum value of ages is used as baseline and n-1 is used instead of n for calculations (Solow 2003). 

# method	Either "StraussSadler" for Strauss & Sadler's (1989) method (default) or "Solow" for Solow's (2003) method. A third option is to use the qbeta function (option "Beta"), since the ratio between the observed timespan and the real timespan for Strauss & Sadler's model is distributed according to a Beta distribution with parameters N and 1 (Wang et al. 2009); this should give the same results as "StraussSadler".


 qage <- function(p=0.5, ages, baseline=NULL, method="StraussSadler") {
 	
 	# If a baseline is not provided, use the minimum of ages as the baseline and subtract one from n.
	
	if(is.null(baseline)) {
		baseline <- min(ages)
		n <- length(ages) -1
		}
		else n <- length(ages)

	# Apply the baseline

	A <- ages - baseline

	if(method=="StraussSadler") {	X <- max(A)*(1-p)^-(1/n) }
	
	else if(method=="Beta") {
		
		Y <- qbeta(p=p, shape1=n, shape2=1)
		
		X <- max(A)/Y}
	
	else if(method=="Solow") {
		
		N <- length(A)
	
		sA <- sort(A)
	
		Z <- sA[N] - sA[N-1]
		
		X <- max(A) + Z*p/(1-p) }
	
	return(X + baseline)
	}


qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5))
qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Beta")
qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Solow")

qage(p=c(0.1, 0.5, 0.9), ages=c(54,51), method="Solow")
qage(p=c(0.1, 0.5, 0.9), ages=c(54,51), method="StraussSadler")

# method="Solow" gives the same answer as method="StraussSadler" when there are only two ages and there is no baseline specified (can be demostrated algebraically).



### Random clade age generator ###

# This function uses the quantile function above to generate random clade ages from the corresponding models.

# Arguments #

# n		Number of numbers to be generated.

# min.ages, max.ages Vector of ages, potentially accounting for age uncertainty. If ages are known exactly, only min.ages is used. If some or all ages have uncertainty, typically upper and lower bounds defined by bracketing geological strata, age maximums and minimums are set for each fóssil.

# max.p		Option for sampling up to a certain probability, usefull for avoiding extremelly large values when using Solow (2003). Usefull for exploratory analysis and histograms but the maximum probability should not be restricted if the sample is used in Monte Carlo simulations.

# ...	Other options passed to qage

rage <- function(n, min.ages, max.ages=NULL, max.p=1, ...) {

	if(is.null(max.ages)) {
		
		rA <- qage(p=runif(n, max=max.p), ages=min.ages, ...)
	}

	else if(length(min.ages)==length(min.ages)) {
	
	rA <- numeric(length=n)
	
	for(i in 1:n) {
	
	A <- runif(length(min.ages), min=min.ages, max=max.ages)
	
	rA[i] <- qage(p=runif(1, max=max.p), ages=A, ...)
		
		}
	}
	
	return(rA)	
}


hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5)), breaks=500, freq=FALSE, xlim=c(50,120))

curve(dexp(x-50, rate=0.07), col='red', lwd=2, add=TRUE)


hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5), max.p=0.95, method="Solow"), breaks=500, freq=FALSE, xlim=c(50,120), border="red")

curve(dexp(x-50, rate=0.04), col='red', lwd=2, add=TRUE) # too heavy on intermediate values

curve(dlnorm(x-50, meanlog=log(50-30), sdlog=pi/sqrt(3)), col='green4', lwd=2, add=TRUE)

# This parameterization was suggested by Norris et al. (manuscript) for the case of using the two oldest fosils (i.e. Solow method): meanlog is set to equal log of last gap; sdlog = pi/sqrt(3).


# With age uncertainty in fossil ages #

hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5), max.ages=c(56, 35, 25, 14, 6)), breaks=500, freq=FALSE, xlim=c(50,120))

curve(dlnorm(x-50, meanlog=2.6, sdlog=0.9), col='red', lwd=2, add=TRUE)

# Since the mode of the empirical density is usually very close to the maximum possible age of the oldest fossil, an alternative parameterization of the log-normal probability density function can facilitate the visual fitting by using the interval range of the oldest fossil to set the mode parameter.

DLnorm <- function(x, mode, sdlog) {
	R <- dlnorm(x, meanlog= sdlog^2+log(mode), sdlog=sdlog)
	cat("meanlog =", sdlog^2+log(mode), "; sdlog =", sdlog=sdlog)
	return(R)
}


hist(sample1, breaks=500, freq=FALSE, xlim=c(50,120))

curve(DLnorm(x-50, sdlog=0.85, mode=56-50), col='red', lwd=3, add=TRUE)


# Age uncertainty and Solow's method

hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5), max.ages=c(56, 35, 25, 14, 6), method="Solow", max.p=0.99), breaks=500, freq=FALSE, xlim=c(50,300))

curve(dlnorm(x-50, meanlog=log(50-30), sdlog=pi/sqrt(3)), col='green4', lwd=2, add=TRUE)

curve(DLnorm(x-50, sdlog=1.1, mode=56-50), col='red', lwd=3, add=TRUE)



# A fossil record with greater age uncertainty

sample2 <- rage(100000, min.ages=c(50, 30, 25, 14, 3.5), max.ages=c(65, 35, 25, 14, 6))

hist(sample2, breaks=1000, freq=FALSE, xlim=c(50,150))

curve(DLnorm(x-50, sdlog=0.6, mode=64-50), col='red', lwd=2, add=TRUE)

curve(dlnorm(x-50, meanlog=2.6, sdlog=0.9), col='red', lwd=2, add=TRUE)


# The greaer the uncertainty in the age of the oldest fósil the greater the difference between the sample and the log-normal function: another advantage of using the sample itself in ML MC instead of a parametric distribution plus Bayesian inference.




### References ###

# Claramunt, S., & J. Cracraft 2015. A new time-tree reveals Earth history's imprint on the evolution of modern birds. ScienceAdvances in press.
# Gingerich, P. D., & M. D. Uhen. 1998. Likelihood estimation of the time of origin of Cetacea and the time of divergence of Cetacea and Artiodactyla. Palaeontologia Electronica 1:45.
# Norris, R. W., C. L. Strope, D. M. McCandlish, A. Stoltzfus (manuscript). Bayesian priors for tree calibration: Evaluating two new approaches based on fossil intervals. bioRxiv, doi: https://doi.org/10.1101/014340
#Solow, A. R. 2003. Estimation of stratigraphic ranges when fossil finds are not randomly distributed. Paleobiology 29(2):181-185.
# Strauss D, & P. M. Sadler 1989. Classical confidence intervals and Bayesian probability estimates for ends of local taxon ranges. Mathematica Geology 21:411-427.
# Wang, S. C. 2010. Principles of statistical inference: likelihood and the Bayesian paradigm. Pp. 1-18 in J. Alroy & G. Hunt (eds.) Quantitative methods in paleobiology. The Paleontological Society Papers 16.
# Wang, S. C., D J. Chudzicki & P. J. Everson 2009. Optimal estimators of the position of a mass extinction when recovery potential is uniform. Paleobiology 35(3):447–459.
# Wang, S. C., & P. J. Everson. 2007. Confidence intervals for pulsed mass extinction events. Paleobiology 33(2):324–336.

