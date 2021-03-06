% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate_rel_change}
\alias{simulate_rel_change}
\title{Simulate concentration changes}
\usage{
simulate_rel_change(n, par1, par2 = NULL, p = 0.5, con = NULL)
}
\arguments{
\item{n}{Number of samples to draw.}

\item{par1}{Parameters (alpha, beta) for first distribution.}

\item{par2}{Optional parameters (alpha, beta) for second distribution.}

\item{p}{Proportion of samples to draw using \code{par1} vs. \code{par2}.
Note, \code{p} is ignored if \code{par2} is not specified.}

\item{con}{Optional vector of lower and upper constraints on generated
samples.}
}
\value{
A vector of relative concentration values.
}
\description{
Simulates relative (fractional) changes in metabolite concentrations using
a mixture of beta distributions.
}
\examples{
out <- simulate_rel_change(10000, c(2, 5), c(0.5, 0.5), 0.7, c(0.1, 1))
hist(out, 20, probability = TRUE,
     main = '', xlab = 'Fractional change in concentration')
}

