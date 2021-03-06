% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/simulation.R
\name{rmix}
\alias{rmix}
\title{Mix random samples}
\usage{
rmix(n, p, rdist1, rdist2)
}
\arguments{
\item{n}{Number of samples to draw.}

\item{p}{Proportion of samples to draw from rdist1 vs rdist2
(true proportion taken as the fraction of random uniform values
 smaller than p).}

\item{rdist1}{User-defined function for generating random samples
(must take \code{n} as only argument).}

\item{rdist2}{User-defined function for generating random samples
(must take \code{n} as only argument).}
}
\value{
A vector of \code{n} values sampled from \code{rdist1} and
        \code{rdis2}.
}
\description{
Generates a random sample by mixing two distributions.
}
\examples{
rdist1 <- function(n) {rnorm(n, -1, 0.5)}
rdist2 <- function(n) {rnorm(n,  1, 0.5)}
out <- rmix(10000, 0.3, rdist1, rdist2)
print(summary(out))
hist(out, 20, xlim = c(-2, 2), probability = TRUE,
     main = '', xlab = 'Random variable value')
}

