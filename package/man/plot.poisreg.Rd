% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.poisreg}
\alias{plot.poisreg}
\title{Plot Trace and Distribution of Regression Parameters}
\usage{
\method{plot}{poisreg}(x, ...)
}
\arguments{
\item{x}{object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The function calls \code{\link[coda]{plot.mcmc}} on the matrix of sampled regression coefficients, and returns the trace of the sampled outputs and a density estimate for each variable in the chain.
}
\description{
Plot Trace and Distribution of Regression Parameters
}
\seealso{
\code{\link{sample_bpr}}, \code{\link[coda]{plot.mcmc}}
}
