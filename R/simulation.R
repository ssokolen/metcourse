# Functions for simulating metabolic time-courses.

#=========================================================================>
# Helper functions

#------------------------------------------------------------------------
#' Constrains Random Samples 
#'
#' Resamples supplied distribution to ensure generated values fall
#' within given constraints.
#
#' @param n Number of samples to draw.
#' @param rdist User-defined function for generating random samples 
#'              (must take \code{n} as only argument).
#' @param con Vector of lower and upper constraints on generated samples.
#' @param max.iter Maximum number of iterations before aborting.
#
#' @return A vector of \code{n} values sampled from \code{rdist} such that
#'         the minimum sampled value is greater than or equal to \code{con[1]}
#'         and the the maximum is smaller than or equal to \code{con[2]}.
#'
#' @examples
#' # Using default distribution parameters
#' out <- rcon(10000, rnorm, c(-1, 1))
#' print(summary(out))
#' hist(out, 20, xlim=c(-2, 2), probability=TRUE, 
#'      main='', xlab='Random variable value')
#'
#' # Modifying mean and standard deviation
#' out <- rcon(10000, function(n) {rnorm(n, 0.5, 0.5)}, c(-1, 1))
#' print(summary(out))
#' hist(out, 20, xlim = c(-2, 2), probability = TRUE, 
#'      main = '', xlab = 'Random variable value')
#' @export
rcon <- function(n, rdist, con, max.iter=100) {

  out <- numeric(n)
  ind <- 1
  n.c <- 0

  for (i in 1:max.iter) {
    
    s <- rdist(n)
    
    logic <- s > con[1] & s < con[2]
    n.s <- min(sum(logic), n-n.c)

    out[ind:(ind + n.s-1)] <- s[logic][1:n.s]
    ind <- ind + n.s
    n.c <- n.c + n.s

    if (n.c == n) break
  }

  if (n.c < n) {
   msg <- 'Maximum number of iterations exceeded.'
   stop(msg) 
  }

  return(out)
}

#------------------------------------------------------------------------
#' Mixes Random Samples 
#'
#' Generates a random sample by mixing two distributions.
#
#' @param n Number of samples to draw.
#' @param p Proportion of samples to draw from rdist1 vs rdist2 
#'          (true proportion taken as the fraction of random uniform values
#'           smaller than p).
#' @param rdist1 User-defined function for generating random samples 
#'              (must take \code{n} as only argument).
#' @param rdist2 User-defined function for generating random samples 
#'              (must take \code{n} as only argument).
#
#' @return A vector of \code{n} values sampled from \code{rdist1} and
#'         \code{rdis2}.
#'
#' @examples
#' rdist1 <- function(n) {rnorm(n, -1, 0.5)}
#' rdist2 <- function(n) {rnorm(n,  1, 0.5)}
#' out <- rmix(10000, 0.3, rdist1, rdist2)
#' print(summary(out))
#' hist(out, 20, xlim = c(-2, 2), probability = TRUE, 
#'      main = '', xlab = 'Random variable value')
#' @export
rmix <- function(n, p, rdist1, rdist2) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  n1 <- ceiling(n*p)
  n2 <- n - n1

  out <- c(rdist1(n1), rdist2(n2))
  out <- sample(out)
  
  return(sample(out, n))
}

#------------------------------------------------------------------------
#' Scales Random Samples
#'
#' Stretches a vector of values between specified constraints.
#
#' @param x Vector of values to stretch.
#' @param con Vector in the form of c(min, max) specifying new constraints 
#
#' @return A vector of values rescaled to have a minimum of \code{con[1]} and
#'         a maximum of \code{con[2]}.
#'
#' @examples
#' out <- stretch(rnorm(10000), c(-2, 1))
#' print(summary(out))
#' hist(out, 20, xlim=c(-2, 2), probability=TRUE, 
#'      main = '', xlab = 'Random variable value')
#' @export
#-------------------------------------------------------------------------
stretch <- function(x, con) {
  x <- (x - min(x))/(max(x) - min(x)) * (con[2] - con[1]) + con[1]

  return(x)
}

#=========================================================================>
# Parameter simulation

#-------------------------------------------------------------------------
#' Simulate Time-course Characteristic
#'
#' Underlying function for simulating maximum concentration, relative
#' changes in concentration and relative standard deviations. Each simulation
#' can be performed with one or a mixture of two underlying distributions.
#'
#' @param n Number of samples to draw.
#' @param dist1 Distribution function for the first distribution e.g. rnorm
#' @param par1 A list of parameters for the first distribution (to be passed 
#'        into dist1 using do.call).
#' @param dist2 Optional distribution function for the second distribution.
#'              Note, \code{par2} has to be specified if \code{dist2} is 
#'              specified.
#' @param par2 An optional list of parameters for the second distribution. 
#'              Note, \code{dist2} has to be specified if \code{par2} is 
#'              specified.
#' @param p Proportion of samples to draw using \code{dist1} vs. \code{dist2}.
#'          Note, \code{p} is ignored if \code{dist2} or \code{par2} 
#'          is not specified.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples.
#
#' @return A vector of randomly generated values.
simulate_characteristic <- function(n, dist1, par1, dist2 = NULL, par2 = NULL, 
                                    p = 0.5, con = NULL) {

  if (xor(is.null(dist2), is.null(par2))) {
    msg <- 'Either both or none of dist2 and par2 must be specified'
    stop(msg)
  }

  if (is.null(dist2)) {
    rtotal <- function(n) do.call(dist1, c(list(n=n), par1))
  } else {
    if ((p < 0) | (p > 1)) {
      msg <- 'p must be between 0 and 1'
      stop(msg)
    } 

    rdist1 <- function(n) do.call(dist1, c(list(n=n), par1))
    rdist2 <- function(n) do.call(dist2, c(list(n=n), par2))
    rtotal <- function(n) rmix(n, p, rdist1, rdist2)
  }

  if (is.null(con)) {
    out <- rtotal(n)
  } else {
    out <- rcon(n, rtotal, con)
  }

  return(out)
}

#-------------------------------------------------------------------------
#' Simulates Maximum Concentrations
#'
#' Simulates maximum metabolite concentrations using a mixture of 2 normal 
#' distributions of concentration logarithms.
#'
#' @param n Number of samples to draw.
#' @param par1 Parameters (mean, sd) in concentration units for first
#'             distribution. These values are converted to logarithms for 
#'             sampling.
#' @param par2 Optional parameters (mean, sd) in concentration units for second
#'             distribution. These values are converted to logarithms for 
#'             sampling.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#'          Note, \code{p} is ignored if \code{par2} is not specified.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples.
#
#' @return A vector of concentration values.
#'
#' @examples
#' # Generating concentrations
#' out <- simulate_max(10000, c(7, 2), c(0.5, 2), 0.3)
#'
#' # Formatting output on logarithmic scale
#' out <- log10(out)
#' labels <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50)
#' hist(out, 20, axes = FALSE, probability = TRUE, 
#'      main = '', xlab = 'Maximum metabolite concentration')
#' axis(side = 2)
#' axis(at = log10(labels), labels = labels, side = 1)
#' @export
simulate_max <- function(n, par1, par2 = NULL, p = 0.5, con = NULL) {

  if (is.null(par2)) {
    if (par1[2] <= 1) {
      msg <- 'Log transformation results in negative or zero standard devation'
      stop(msg)
    }

    out <- simulate_characteristic(n, rnorm, as.list(log10(par1)), con=con)

  } else {
    if ((par1[2] <= 1) | par2[2] <= 1) {
      msg <- 'Log transformation results in negative or zero standard devation'
      stop(msg)
    }

  if (!is.null(con)) {
    con <- log10(con)
  }

  out <- simulate_characteristic(n, rnorm, as.list(log10(par1)), 
                                    rnorm, as.list(log10(par2)),                                       
                                    p, con)
  }

  return(10**out)
}

#-------------------------------------------------------------------------
#' Simulate Concentration Changes 
#'
#' Simulates relative (fractional) changes in metabolite concentrations using
#' a mixture of beta distributions.
#'
#' @param n Number of samples to draw.
#' @param par1 Parameters (alpha, beta) for first distribution.
#' @param par2 Optional parameters (alpha, beta) for second distribution.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#'          Note, \code{p} is ignored if \code{par2} is not specified.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples.
#
#' @return A vector of relative concentration values.
#'
#' @examples
#' out <- simulate_rel_change(10000, c(2, 5), c(0.5, 0.5), 0.7, c(0.1, 1))
#' hist(out, 20, probability = TRUE, 
#'      main = '', xlab = 'Fractional change in concentration')
#' @export
simulate_rel_change <- function(n, par1, par2 = NULL, p = 0.5, con = NULL) {

  if (is.null(par2)) {
    out <- simulate_characteristic(n, rbeta, as.list(par1), con=con)

  } else {
    out <- simulate_characteristic(n, rbeta, as.list(par1), 
                                      rbeta, as.list(par2), 
                                      p, con)
  }

  return(out)
}

#-------------------------------------------------------------------------
#' Simulate Measurement Variability 
#'
#' Simulates relative standard deviations (coefficients of variation) using
#' a mixture of 2 normal distributions (one subpopulation of less variable
#' measurements and on subpopulation of more variable ones).
#'
#' @param n Number of samples to draw.
#' @param par1 Parameters (mean, sd) for first distribution.
#' @param par2 Optional parameters (mean, sd) for second distribution.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#'          Note, \code{p} is ignored if \code{par2} is not specified.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples. By default, relative standard deviations are limited
#'            to between 0 and 1.
#
#' @return A vector of relative standard deviation values.
#'
#' @examples
#' out <- simulate_rel_sd(10000, c(.04, .02), c(0.11, 0.02), 0.7, c(0, 0.20))
#' hist(out, 20, probability = TRUE, 
#'      main = '', xlab = 'Relative standard deviation')
#' @export
simulate_rel_sd <- function(n, par1, par2 = NULL, p = 0.05, con = c(0, 1)) {

  if (is.null(par2)) {
    out <- simulate_characteristic(n, rnorm, as.list(par1), con=con)

  } else {
    out <- simulate_characteristic(n, rnorm, as.list(par1), 
                                      rnorm, as.list(par2), 
                                      p, con)
  }

  return(out)
}

#=========================================================================>
# Trend simulation

#-------------------------------------------------------------------------
#' Simulate Trends
#'
#' Underlying function for simulating different types of time-course trends.
#' Each trend supports parameter generation from one or two sets of parameter
#' limits.
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param f_trend Function that taxes x values and a vector of parameters
#'                and generates a corresponding set of y values.
#' @param par1 Parameter vector of lower and upper bounds for first set of 
#'        trends.
#' @param par2 Optional parameter vector of lower and upper bounds for second 
#'        set of trends.
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#'          Note, \code{p} is ignored if \code{par2} is not specified.
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
simulate_trend <- function(n, n.samples, f_trend, par1, par2 = NULL, p = 0.5) {

  if (is.null(par2)) {
    n1 <- n
    n2 <- 0
  } else {
    if ((p < 0) | (p > 1)) {
      msg <- 'p must be between 0 and 1'
      stop(msg)
    } 

    n1 <- floor(n*p)
    n2 <- n - n1
  }

  # Initializing parameter matrix
  n.par <- length(par1)
  par <- matrix(NA, nrow=n, ncol=n.par/2)

  # Generating trends from first set of parameters
  if (n1 != 0) {
    for (i in seq(1, n.par, by=2)) {
      par[1:n1, (i+1)/2] <- runif(n1, par1[i], par1[i+1])
    }
  }

  # Generating trends from second set of parameters
  if (n2 != 0) {
    for (i in seq(1, n.par, by=2)) {
      par[(n1+1):n, (i+1)/2] <- runif(n2, par2[i], par2[i+1])
    }
  }

  # Shuffling order
  par <- par[sample(1:n), ]

  # Defining x variables
  x <- seq(0, 1, length.out=n.samples)

  # Forcing 1-row matrix to be a matrix
  if (is.null(dim(par))) {
    dim(par) <- c(1, n.par/2)
  }

  # Applying f_trend on every set of parameters
  out <- apply(par, 1, function(par) f_trend(x, par))

  # Taking the transpose to get samples as columns
  out <- t(out)

  return(as.data.frame(out))
}

#-------------------------------------------------------------------------
#' Simulate Decreasing Trends
#'
#' Simulates decreasing trends by mixing two sigmoid curve families. Sigmoid
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function. Function: y = (1 + exp((x-a)/b))**(-1).
#'
#' 
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max) for first set of 
#'        trends
#' @param par2 Optional parameter vector (a_min, a_max, b_min, b_max) for second 
#'        set of trends
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#'        Note, \code{p} is ignored if \code{par2} is not specified.
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(0.2, 0.6, 0.10, 0.18)
#' par2 <- c(0.6, 0.9, 0.10, 0.18)
#' trends <- simulate_decreasing(1000, 100, par1, par2, 0.05)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_decreasing <- function(n, n.samples, par1, par2 = NULL, p = 0.5) {

  # Defining function
  f_dec <- function(x, par) {
    y <- (1 + exp((x - par[1]) / par[2]))**(-1)
    out <- stretch(y, c(0, 1))
    return(out)
  }

  out <- simulate_trend(n, n.samples, f_dec, par1, par2, p)

  return(out)
}

#-------------------------------------------------------------------------
#' Simulate Increasing Trends
#'
#' Simulates increasing trends by mixing two sigmoid curve families. Sigmoid
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function. Function: y = 1 - (1 + exp((x-a)/b))**(-1).
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max) for first set of 
#'        trends
#' @param par2 Optional parameter vector (a_min, a_max, b_min, b_max) for second 
#'        set of trends
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#'        Note, \code{p} is ignored if \code{par2} is not specified.
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(0.045, 0.055, 0.2, 0.4)
#' par2 <- c(0.945, 0.955, 0.1, 0.3)
#' trends <- simulate_increasing(1000, 100, par1, par2, 0.15)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_increasing <- function(n, n.samples, par1, par2 = NULL, p = 0.5) {

  # Defining function
  f_inc <- function(x, par) {
    y <- 1 - (1 + exp((x - par[1]) / par[2]))**(-1)
    out <- stretch(y, c(0, 1))
    return(out)
  }

  out <- simulate_trend(n, n.samples, f_inc, par1, par2, p)

  return(out)
}

#-------------------------------------------------------------------------
#' Simulate Concave Trends
#'
#' Simulates concave trends by mixing trimmed and scaled beta distributions. 
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function. Function y = x**(a-1)*(1-x)**(b-1), with the
#' domain trimmed such that min(x) = c and max(x) = d (before rescaling 
#' to [0, 1])
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max) for first set of 
#'        trends
#' @param par2 Optional parameter vector (a_min, a_max, b_min, b_max) for second 
#'        set of trends
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#'        Note, \code{p} is ignored if \code{par2} is not specified.
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9)
#' par2 <- c(3.5, 4.5, 2.5, 3.5, 0.2, 0.4, 0.7, 0.8)
#' trends <- simulate_concave(1000, 100, par1, par2, 0.75)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_concave <- function(n, n.samples, par1, par2 = NULL, p = 0.5) {

  # Defining function
  f_conc <- function(x, par) {

    # Sampling at higher density to avoid problems with trimming later
    x.temp <- seq(min(x), max(x), length.out=length(x)*10)
    y <- x.temp**(par[1] - 1) * (1 - x.temp)**(par[2] - 1)

    # Trimming sections of the curve as specified
    logic <- (x.temp > par[3]) & (x.temp < par[4])
    y <- y[logic]
    x.temp <- x.temp[logic]

    # Rescaling to [0, 1]
    y <- stretch(y, c(0, 1))
    x.temp <- stretch(x.temp, c(0, 1))

    f_spline <- splinefun(x.temp, y)

    # Generating y values at initially desired x
    out <- f_spline(x)
    return(out)
  }

  out <- simulate_trend(n, n.samples, f_conc, par1, par2, p)

  return(out)
} 

#-------------------------------------------------------------------------
#' Simulate Linear Trends
#'
#' Simulates linear trends with p fraction decreasing and (1-p) fraction
#' increasing.
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param p Proportion of trends decreasing (vs. increasing).
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' trends <- simulate_linear(1000, 100, 0.75)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_linear <- function(n, n.samples, p) {

  n1 <- floor(n*p)
  n2 <- n - n1

  # Defining function
  x <- seq(0, 1, length.out=n.samples)
  
  y1 <- seq(1, 0, length.out=n.samples)
  y2 <- x

  out1 <- matrix(rep(y1, each=n1), nrow=n1, ncol=n.samples)
  out2 <- matrix(rep(y2, each=n2), nrow=n2, ncol=n.samples)

  out <- rbind(out1, out2)
  out <- out[sample(1:n), ]

  if (is.null(dim(out))) {
    dim(out) <- c(1, n.samples)
  }

  return(as.data.frame(out))
} 

#-------------------------------------------------------------------------
#' Simulate Full Metabolic Time-courses
#'
#' Combines various simulating functions to generate "experiments" composed
#' of multiple trends based on provided parameters. This function is meant 
#' as a template for custom simulations.
#'
#' @param n.trends Number of trends to generate per "experiment" i.e. number
#'                 of metabolites.
#' @param p A vector of fractions representing the proportions of decreasing, 
#'          increasing, and concave trends (in that order). If the fractions
#'          don't sum to one, the remainder of the trends are assumed to be
#'          linear (undetermined). Note that floor() will be applied to convert
#'          final numbers to integers.
#' @param param A list of parameters to be fed into the simulation of various
#'              time-course components. Parameters must be specified for
#'              maximum concentration, relative changes in concentration,
#'              measurement standard deviations, and the trends themselves.
#'              Each component has a qualifier (max, change, sd, trend) that is
#'              appended to the name of the argument used in the function 
#'              e.g. max.par1 is used to set par1 for simulating maximum
#'              concentrations. To specify different parameters for each type
#'              of trend, append trend qualifier e.g. max.par1.decreasing can
#'              be used to set the maximum concentrations of only decreasing
#'              trends, with max.par1 used for all other trends.
#' @param n.samples Number of timepoints within each trend.
#' @param n.experiments Number of "experiments".
#
#' @return A "long" dataframe with the following columns: experiment, 
#'         metabolite, sample, concentration.
#'
#' @examples
#' # Constructing realistic list of parameters
#' param <- list(
#'   # Maximum concentrations are the same for every trend type
#'   p.max = 0.3, 
#'   par1.max = c(7, 2), 
#'   par2.max = c(0.5, 2),  
#'   con.max = c(0, 50), 
#'               
#'   # Global change parameters are near 100% for increasing/concave trends
#'   par1.change = c(5, 0.1), 
#'   con.change = c(0.5, 1),
#' 
#'   # Decreasing trends can have a wide variety of changes
#'   p.change.decreasing = 0.7,  
#'   par1.change.decreasing = c(2, 5), 
#'   par2.change.decreasing = c(0.5, 0.5),
#'   con.change.decreasing = c(0.1, 1),
#' 
#'   # Linear trends are characterized by relatively small changes
#'   par1.change.linear = c(1, 5), 
#'   con.change.linear = c(0, 0.1),
#' 
#'   # Measurement error is the same for every trend type (but no more than 20%)
#'   p.sd = 0.7, 
#'   par1.sd = c(0.04, 0.02), 
#'   par2.sd = c(0.11, 0.02),
#'   con.sd = c(0, 0.20),
#' 
#'   # Decreasing trend specification
#'   p.trend.decreasing = 0.05,
#'   par1.trend.decreasing = c(0.2, 0.6, 0.10, 0.18),
#'   par2.trend.decreasing = c(0.6, 0.9, 0.10, 0.18),
#' 
#'   # Increasing trend specification
#'   p.trend.increasing = 0.15,
#'   par1.trend.increasing = c(0.045, 0.055, 0.2, 0.4),
#'   par2.trend.increasing = c(0.945, 0.955, 0.1, 0.3),
#' 
#'   # Concave trend specification
#'   par1.trend.concave = c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9),
#' 
#'   # Linear trends are equaly split between increasing and decreasing
#'   p.trend.linear = 0.5
#' )
#'
#' # Generating trends
#' timecourse <- simulate_timecourse(10, c(0.3, 0.3, 0.3), param)
#'
#' # Plotting
#' par(mfrow = c(5, 2), oma = c(5, 4, 1, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)
#' 
#' for (metabolite in unique(timecourse$metabolite)) {
#'   logic <- timecourse$metabolite == metabolite  
#'   plot(timecourse$sample[logic], timecourse$concentration[logic],
#'        xlab='', ylab='')
#' }
#'
#' title(xlab = 'Sample', ylab = 'Concentration', outer = TRUE, line = 3)
#' @export
simulate_timecourse <- function(n.trends, p, param,
                                n.samples = 10, n.experiments = 1) {

  if (any((p < 0) | (p > 1))) {
   msg <- 'All elements of p must be between 0 and 1'
   stop(msg)
  } 

  if (sum(p) > 1) {
   msg <- 'Elements of p must sum to less than or equal to 1'
   stop(msg)
  } 

  # Converting fractions to numbers
  n.total <- n.trends
  n.trends <- floor(p*n.trends)
  n.trends <- c(n.trends, n.total - sum(n.trends))
  names(n.trends) <- c('decreasing', 'increasing', 'concave', 'linear')
  n.trends <- as.list(n.trends)

  ### Parsing parameters

  # Checking input
  valid.parameters <- c('p', 'par1', 'par2', 'con')
  valid.simulations <- c('max', 'change', 'sd', 'trend')
  valid.trends <- c('decreasing', 'increasing', 'concave', 'linear')

  general.param <- apply(expand.grid(valid.parameters, 
                                     valid.simulations),
                         1, paste, collapse='.')
  specific.param <- apply(expand.grid(valid.parameters, 
                                      valid.simulations,
                                      valid.trends),
                         1, paste, collapse='.')

  param.names <- names(param)
  invalid <- param.names[! param.names %in% c(general.param, specific.param)]

  if (length(invalid) > 0 ) {
    msg <- sprintf('The following parameters are not valid: %s',
                   paste(invalid, collapse=', '))
    stop(msg)
  }
  
  # Initializing hierarchy
  entries <- list('max'=list(), 'change'=list(), 'sd'=list(), 'trend'=list())
  new.param <- list('decreasing'=entries, 'increasing'=entries,
                    'concave'=entries, 'linear'=entries)

  # Initializing all parameters to NULL  
  for (entry in specific.param) {
    tags <- strsplit(entry, '\\.')[[1]]
    new.param[[tags[3]]][[tags[2]]][[tags[1]]] <- NULL
  }

  # Filling in general parameters 
  for (entry in param.names[param.names %in% general.param]) {
    tags <- strsplit(entry, '\\.')[[1]]
    
    for (trend in valid.trends) {
      new.param[[trend]][[tags[2]]][[tags[1]]] <- param[[entry]]
    }
  }

  # Filling in specific parameters
  for (entry in param.names[param.names %in% specific.param]) {
    tags <- strsplit(entry, '\\.')[[1]]
    new.param[[tags[3]]][[tags[2]]][[tags[1]]] <- param[[entry]]
  }

  ### Generating timecourses
  simulate_trend <- list('decreasing' = simulate_decreasing, 
                         'increasing' = simulate_increasing,
                         'concave' = simulate_concave,
                         'linear' = simulate_linear)

  # Initializing global data frame
  combined <- c()

  for (trend in valid.trends) {

    # Generating experiment and metabolite tags
    experiments <- rep(1:n.experiments, each=n.trends[[trend]])
    metabolites <- paste(rep(1:n.trends[[trend]], n.experiments), trend, sep='')

    # Calculating simulation number
    n <- n.trends[[trend]] * n.experiments

    # Generating trends
    simulation <- 'trend'
    par1 <- new.param[[trend]][[simulation]][['par1']]
    par2 <- new.param[[trend]][[simulation]][['par2']]
    p <- new.param[[trend]][[simulation]][['p']]
    
    if (trend == 'linear') {
      trends <- simulate_trend[[trend]](n, n.samples, p)
    } else {
      trends <- simulate_trend[[trend]](n, n.samples, par1, par2, p)
    }

    sample.names <- paste('s', 1:n.samples, sep = '')
    colnames(trends) <- sample.names

    trends$experiment <- experiments
    trends$metabolite <- metabolites

    # Initializing single parameter dataframe
    parameters <- data.frame(experiment=experiments, metabolite=metabolites,
                             stringsAsFactors=FALSE)

    # Maximum concentrations
    simulation <- 'max'
    p <- new.param[[trend]][[simulation]][['p']]
    par1 <- new.param[[trend]][[simulation]][['par1']]
    par2 <- new.param[[trend]][[simulation]][['par2']]
    con <- new.param[[trend]][[simulation]][['con']]
    parameters$max <- simulate_max(n, par1, par2, p, con)

    # Minimum concentrations (calculated from percent change)
    simulation <- 'change'
    p <- new.param[[trend]][[simulation]][['p']]
    par1 <- new.param[[trend]][[simulation]][['par1']]
    par2 <- new.param[[trend]][[simulation]][['par2']]
    con <- new.param[[trend]][[simulation]][['con']]
    changes <- simulate_rel_change(n, par1, par2, p, con)
    parameters$min <- parameters$max * (1 - changes)

    # Measurement standard deviations
    simulation <- 'sd'
    p <- new.param[[trend]][[simulation]][['p']]
    par1 <- new.param[[trend]][[simulation]][['par1']]
    par2 <- new.param[[trend]][[simulation]][['par2']]
    con <- new.param[[trend]][[simulation]][['con']]
    parameters$sd <- simulate_rel_sd(n, par1, par2, p, con)

    # Combining
    trends <- left_join(trends, parameters, by = c('experiment', 'metabolite'))

    # Converting to long format
    trends <- gather_(trends, 'sample', 'conc', sample.names)

    # Applying parameters
    trends <- trends %>%
                group_by(experiment, metabolite) %>%
                mutate(conc = conc * (max - min) + min,
                       conc = conc + rnorm(n(), 0, mean(conc) * sd)) %>%
                select(experiment, metabolite, sample, conc) %>%
                rename(concentration=conc) %>%
                ungroup()

    combined <- rbind(combined, trends)
  }

  # Randomizing metabolite names
  met.names <- sample(1:n.total)
  names(met.names) <- unique(combined$metabolite)

  combined$metabolite <- as.numeric(met.names[combined$metabolite])
  combined$sample <- as.numeric(gsub('s', '', combined$sample))
  combined$experiment <- as.numeric(combined$experiment)

  combined <- arrange(combined, experiment, metabolite, sample)

  return(combined)
}
