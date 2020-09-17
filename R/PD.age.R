
##' @title Probability densities of clade age estimates derived from descendant
##' fossil occurrences
##' @description \code{PD.age} assumes that fossil ages are known without error.
##' In contrast, \code{PD.age.unc} incorporates fossil age uncertainty; in the
##' latter, the minimum and maximum ages are treated as bounds on a uniform
##' distribution. \code{PD.age.plot} plots the distribution of interest.
##' @param x A vector of fossil ages (\code{PD.age}) or a matrix
##' (\code{PD.age.unc}) of fossil ages in which the first column is the maximum
##' age (or lower stratigraphic bound) and the second column is the minimum
##' fossil age (upper stratigraphic bound) for each fossil.
##' @param baseline Youngest bound of the distribution that could be the
##' present (=0), or some other stratum of reference. If a baseline is not
##' provided, the minimum value of ages is used as baseline and n-1 is used for
##' the sample size instead of n for calculations (Solow, 2003).
##' @param p.max The cumulative probability value to set a limit on the range of
##' values for which the probability is calculated
##' @param reps The number of pseudoreplicates used to account for fossil age
##' uncertainty
##' @param breaks The number of breaks or time intervals used to summarize the
##' pseudoreplicates
##' @param line.col The colour of the probability distribution curve. Only used
##' with \code{PD.age.plot}
##' @param ... Further arguments passed to \code{PD.age.plot}
##' @return In the case of \code{PD.age} and \code{PD.age.unc}, a dataframe with
##' two columns: \code{age} (the geological age under consideration), and
##' \code{p} (the associated probability for the clade origin).
##' @importFrom grDevices gray
##' @importFrom graphics abline axis mtext rect box lines par plot
##' @examples
##' \dontrun{
##'   # No fossil age uncertainty
##'   PD.age(c(54,31, 25, 14, 5));
##'   PD.age(c(54,31, 25, 14, 5), baseline=0);
##'   
##'   # With fossil age uncertainty
##'   PD.age.unc(cbind(c(56, 35, 25, 14, 6), c(50, 30, 25, 14, 3.5)), baseline=0);
##'   
##'   # Plot and fit of probability density function
##'   PD.age.plot(c(54,31, 25, 14, 5), p.max=0.99);
##'   curve(dexp(x-54, rate=0.062), col='black', lwd=3, add=TRUE);
##'   
##'   # Same, but use lognormal distribution
##'   PD.age.plot(cbind(c(56, 35, 25, 14, 6), c(50, 30, 25, 14, 3.5)));
##'   curve(dlnorm(x-50, meanlog=2.7, sdlog=1)*1.3, col='black', lwd=3, add=TRUE);
##'   }
##' @importFrom Rdpack reprompt
##' @references
##' \insertRef{Claramunt2015}{cladeage}
##' 
##' \insertRef{Solow2003}{cladeage}
##' 
##' \insertRef{Wang2007}{cladeage}
##' 
##' \insertRef{Wang2010}{cladeage}
##' @export
PD.age <- function(x, baseline=NULL, p.max=0.99) {
  
  # If a baseline is not provided, use the minimum of x as the baseline and
  # subtract one from n.
  if (is.null(baseline)) {
    baseline <- min(x);
    n <- length(x) -1;
  } else {
    n <- length(x);
  }
  
  # Apply the baseline
  X <- x - baseline;
  
  # Calculate the maximum value to be evaluated, set as the value at the p.max
  # cumulative probability. See Solow (2003) equation 4. 
  MAX <- max(X)*(1-p.max)^-(1/n);
  
  # Generate a sequence of ages for which the probability will be calculated
  Y <- max(X):MAX;
  
  # Calculate the likelihood of those ages Yi given the observed fossil ages,
  # which is numerically equal to the joint probability of observing X given
  # that the age of the clade is Y (Wang 2010):
  Lik <- 1/Y^n;
  
  # Calculate the area under the curve. Because time intervals are 1 and it is a
  # monotonically decreasing curve
  area <- sum(Lik) - Lik[1]/2;
  
  # Rescale Likelihood so area = p.max
  P <- Lik*p.max/area;
  pd <- as.data.frame(cbind(x=Y+baseline, P));
  colnames(pd) <- c("age", "p");
  
  return (pd);
}

##' @export
##' @rdname PD.age
PD.age.unc <- function(x, baseline=NULL, p.max=0.99, reps=1000, breaks=100) {
  
  N <- nrow(x);
  
  lower <- x[,1];
  upper <- x[,2];
  
  simages <- matrix(nrow=N, ncol=reps)
  for (i in 1:N) {
    simages[i,] <- runif(reps, min= upper[i], max= lower[i]);
  }
  
  slist <- list();
  for (i in 1:reps) {
      slist[[i]] <- simages[,i];
  }
  
  UPDFs <- lapply(slist, PD.age, baseline=baseline, p.max=p.max)
  
  # Put all values in two vectors
  UPDFsx <- numeric();
  UPDFsPD <- numeric();
  for (i in 1:length(UPDFs)) {
    UPDFsx <- append(UPDFsx, UPDFs[[i]][,'age']);
    UPDFsPD <- append(UPDFsPD, UPDFs[[i]][,'p']);
  }
  
  Breaks <- seq(from=max(upper), to=max(UPDFsx), length.out=breaks);
  Density <- numeric();
  
  for (i in 2:length(Breaks)-1) {
    Density[i] <- sum(UPDFsPD[UPDFsx >= Breaks[i] & UPDFsx < Breaks[i+1]]);
  }
  
  Density <- Density/reps;
  Results <- as.data.frame(cbind(age=Breaks[-length(Breaks)], p=Density));
  
  return (Results);
}

##' @export
##' @rdname PD.age
PD.age.plot <- function(x, baseline=NULL, p.max=0.99, reps=1000, breaks=100,
                        line.col="red", ...) {
  
  x <- as.matrix(x);
  graphics::par(las=1, lend=1, ljoin=1, yaxs='i');
  
  if (dim(x)[2] == 1) {
    dens <- PD.age(x, baseline=baseline, p.max=p.max);
    
    graphics::plot(dens, type='n', ylim=c(0, max(dens$p)*1.1), ann=FALSE,
                   axes=FALSE, ...);
    graphics::box();
    graphics::abline(v=max(x), lwd=2, col=grDevices::gray(0.8));
    graphics::lines(dens, type='l', lwd=5, col=line.col);
    graphics::axis(1);
    graphics::mtext("time", 1, line=2.5);
    graphics::mtext("P", 2, las=1, line=2);
  } else if (dim(x)[2] == 2) {
    dens <- PD.age.unc(x, baseline=baseline, p.max=p.max, reps=reps,
                       breaks=breaks);
    
    graphics::plot(dens, type='n', ylim=c(0,max(dens$p)*1.1), axes=FALSE,
                   ann=FALSE, ...);
    graphics::box();
    graphics::rect(xleft=max(x[,2]), ybottom=0, xright=max(x[,1]), ytop=1,
                   col=grDevices::gray(0.8), border=grDevices::gray(0.8), lwd=2);
    graphics::lines(dens, type="s", lwd=5, col=line.col);
    graphics::axis(1);
    graphics::mtext("time", 1, line=2.5);
    graphics::mtext("P", 2, las=1, line=2);
  } else {
    cat("Error: check the format of input data x");
  }
}
