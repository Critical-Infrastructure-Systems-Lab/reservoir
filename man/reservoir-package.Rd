\name{reservoir-package}
\alias{reservoir-package}
\alias{reservoir}
\docType{package}
\title{An R package for the analysis, design, and operation of water supply reservoirs}
\description{Measure single reservoir performance using resilience, reliability,
    and vulnerability metrics; compute storage-yield-reliability relationships; 
    determine no-fail Rippl storage with sequent peak analysis; optimize release
    decisions using determinisitc and stochastic dynamic programming; evaluate
    inflow characteristics.}
\details{
\tabular{ll}{
Package: \tab reservoir\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2015-06-17\cr
License: \tab GPL (>= 2)\cr
}
}
\seealso{ Functions: \code{\link{Rippl}}, \code{\link{storage}}, \code{\link{yield}}, \code{\link{rrv}}, \code{\link{dp}}, \code{\link{sdp}}, \code{\link{Hurst}}.}
\author{Sean Turner <swd.turner@gmail.com> and Stefano Galelli}
\keyword{ package }
\examples{
# 1. Express the distribution of Rippl storage for a known inflow process
# a) Assume the inflow process follows a lognormal distribution
# (meanlog = 0, sdlog = 1):
x <- rlnorm(1200)
#
# b) Convert to a monthly time series object beginning Jan 1900
# (ts objects are not required for most functions in "reservoir")
x <- ts(x, start = c(1900, 1), frequency = 12)
#
# c) Begin reservoir analysis... e.g., compute the Rippl storage
x_Rippl <- Rippl(x, R = rep(mean(x) * 0.7, length = length(x)))
no_fail_storage <- x_Rippl$Rippl_storage
#
# d) Resample x and loop the procedure multiple times to get the
# distribution of no-failure storage for the inflow process and 
# constant release (R) equal to 70 percent of the mean inflow.
no_fail_storage <- vector("numeric", 500)
for (i in 1:length(no_fail_storage)){
x <- ts(rlnorm(1200), start = c(1900, 1), frequency = 12)
no_fail_storage[i] <- Rippl(x, R = rep(mean(x) * 0.7,
length = length(x)),plot=FALSE)$Rippl_storage
}
hist(no_fail_storage)
}

