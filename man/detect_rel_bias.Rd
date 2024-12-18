% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detection.R
\name{detect_rel_bias}
\alias{detect_rel_bias}
\title{Detect systematic bias in metabolomic time-course data}
\usage{
detect_rel_bias(time, concentration, metabolite, min.deviation = NULL, ...)
}
\arguments{
\item{time}{A vector of times or sample numbers for metabolite time-courses.
Note, there should be a time value for every \code{concentration}
for every \code{metabolite} e.g. time = c(1, 2, 3, 1, 2, 3),
concentration = c(20, 15, 10, 3, 6, 9),
metabolite = c('glc', 'glc', 'glc', 'lac', 'lac', 'lac').}

\item{concentration}{A vector of metabolite concentrations.}

\item{metabolite}{A vector of metabolite names that correspond to \code{time}
and \code{concentration}.}

\item{min.deviation}{Smallest median relative deviation to identify as a
bias. If not supplied, it will be estimated as 50% of the underlying noise in the data.}

\item{...}{Arguments to be passed into \code{fit_b_spline}, (such as degree and number of knots).}
}
\value{
A dataframe corresponding to the time points with and without systematic error.
}
\description{
Identifies systematic deviations in metabolic data using a B-spline fit.
Timepoints that have a median relative deviation above a threshold value
are assumed to be influenced by a measurement or methodological bias as
compared to the overall trends metabolite concentrations.
}
\examples{

# Using previously simulated data 40 metabolic trends with 10 time points
data(timecourse)

# Adding an error of 5\% at sample 4
logic <- timecourse$sample == 4
timecourse$concentration[logic] <- timecourse$concentration[logic] * 1.05

# Estimating
error <- correct_rel_bias(timecourse$time,
                                         timecourse$concentration,
                                         timecourse$metabolite)


}
