
##' @title Random clade age generator
##' @description This function uses the quantile function above to generate
##' random clade ages from the corresponding models.
##' @param n Number of samples to be frawn
##' @param min.ages A vector of minimum ages, potentially accounting for age
##' uncertainty.
##' @param max.ages A vector of minimum ages, potentially accounting for age
##' uncertainty.
##' @param max.p Option for sampling up to a certain probability, useful for
##' avoiding extremely large values when using Solow (2003). Useful for
##' exploratory analysis and histograms but the maximum probability should not
##' be restricted if the sample is used in Monte Carlo simulations.
##' @param ... Other options passed to \code{qage}
##' @details If ages are known exactly, only min.ages is used. If some or all
##' ages have uncertainty, typically upper and lower bounds defined by
##' bracketing geological strata, age maximums and minimums are set for each
##' fossil.
##' @return A numeric vector of length \code{n} representing simulated clade
##' ages
##' @importFrom stats runif
##' @examples
##' \dontrun{
##'   # Example using the default method of Strauss & Sadler (1989)
##'   hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5)), breaks=500,
##'     freq=FALSE, xlim=c(50,120), xlab="Clade Age",
##'     main="Strauss & Sadler Example\nn=1000, min.ages=c(50, 30, 25, 14, 3.5), max.p=0.95");
##'   curve(dexp(x-50, rate=0.07), col='blue', lwd=2, add=TRUE);
##'   legend("right", legend="dexp", col="blue", lty=1, lwd=2, box.lty=0);
##'   
##'   # Example using the method of Solow (2003)
##'   hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5), max.p=0.95,
##'     method="Solow"), breaks=500, freq=FALSE, xlim=c(50,120), xlab="Clade Age",
##'     main="Solow Example\nn=1000, min.ages=c(50, 30, 25, 14, 3.5), max.p=0.95");
##'   curve(dexp(x-50, rate=0.04), col='blue', lwd=2, add=TRUE); # heavy on intermediate values
##'   curve(dlnorm(x-50, meanlog=log(50-30), sdlog=pi/sqrt(3)), col='orange', lwd=2, add=TRUE);
##'   legend("right", legend=c("dexp", "dlnorm"), col=c("blue", "orange"), lty=1, lwd=2, box.lty=0);
##'   
##'   # Strauss & Sadler (1989) method incorporating fossil age uncertainty
##'   hist(rage(n=10000, min.ages=c(50, 30, 25, 14, 3.5),
##'     max.ages=c(56, 35, 25, 14, 6)), breaks=500, freq=FALSE, xlim=c(50,120),
##'     xlab="Age", main="Solow Example (with fossil age uncertainty)\nn=1000,
##'     min.ages=c(50, 30, 25, 14, 3.5), max.ages=c(56, 35, 25, 14, 6))\nmax.p=0.95");
##'   curve(dlnorm(x-50, meanlog=2.6, sdlog=0.9), col='blue', lwd=2, add=TRUE);
##'   legend("right", legend="dlnorm", col="blue", lty=1, lwd=2, box.lty=0);
##'   }
##' @export
rage <- function(n, min.ages, max.ages=NULL, max.p=1, ...) {
  if (is.null(max.ages)) {
    rA <- qage(p=stats::runif(n, max=max.p), ages=min.ages, ...);
  } else if (length(min.ages)==length(min.ages)) {
    rA <- numeric(length=n);
    for (i in 1:n) {
      A <- stats::runif(length(min.ages), min=min.ages, max=max.ages);
      rA[i] <- qage(p=runif(1, max=max.p), ages=A, ...);
    }
  }
  return (rA);
}
