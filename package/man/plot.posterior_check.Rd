% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ppc.R
\name{plot.posterior_check}
\alias{plot.posterior_check}
\title{Graphical Posterior Predictive Checks}
\usage{
\method{plot}{posterior_check}(x, ...)
}
\arguments{
\item{x}{object of class "\code{posterior_check}" (usually, the result of a call to \code{\link{posterior_predictive}}).}

\item{...}{other parameters to be passed through to plotting functions. See Details.}
}
\value{
The function outputs (at least) three plots for graphical posterior predictive check.\cr
The first plot compares the empirical cumulative distribution function (ECDF) with the cumulative distribution function obtained 
with samples from the posterior predictive distribution (median and point-wise 95\% credible bands). \cr
The second plot compares the distribution of the observed sample with the predictive distribution obtained using the maximum a posteriori (MAP) 
estimates of the regression parameters.\cr
The third plot compares the predictive distribution of a statistic (default is the mean) with the observed value of the same statistics,
 displayed with a red line.
}
\description{
This function is a method for class \code{posterior_check}. Plot diagnostic statistics for graphical posterior predictive checks.
}
\details{
It is possible to generate additional plots that compare the posterior predictive distribution of a statistics with the observed value. 
This is done through the parameter \code{stats}: it is a list with elements the function names of the statistics one wants to compare. 
Default is \code{stats = list("mean")}, other possible values are, e.g., "median", "sd", "max" etc.
}
\examples{
library(MASS) # load the data set
head(epil)

fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
                   iter = 1000)
\donttest{ plot(posterior_predictive(fit), stats = c("mean", "sd", "max"))   }
# plots for posterior predictive check


}
\seealso{
\code{\link{posterior_predictive}}
}
