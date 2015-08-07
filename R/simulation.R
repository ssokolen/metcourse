# Functions for simulating metabolic time-courses.

#=========================================================================>
# Helper functions

#------------------------------------------------------------------------
#' Constraining random samples 
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
#' Mixing random samples 
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

  n1 <- floor(n*p)
  n2 <- n - n1

  out <- c(rdist1(n1), rdist2(n2))
  out <- sample(out)
  
  return(sample(out, n))
}

#------------------------------------------------------------------------
#' Scaling random samples
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
#' Simulating maximum concentrations
#'
#' Simulates maximum metabolite concentrations using a mixture of 2 normal 
#' distributions of concentration logarithms.
#'
#' @param n Number of samples to draw.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#' @param par1 Parameters (mean, sd) in concentration units for first
#'             distribution. These values are converted to logarithms for 
#'             sampling.
#' @param par2 Parameters (mean, sd) in concentration units for second
#'             distribution. These values are converted to logarithms for 
#'             sampling.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples.
#
#' @return A vector of concentration values.
#'
#' @examples
#' # Generating concentrations
#' out <- simulate_max(10000, 0.3, c(7, 2), c(0.5, 2))
#'
#' # Formatting output on logarithmic scale
#' out <- log10(out)
#' labels <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50)
#' hist(out, 20, axes = FALSE, probability = TRUE, 
#'      main = '', xlab = 'Maximum metabolite concentration')
#' axis(side = 2)
#' axis(at = log10(labels), labels = labels, side = 1)
#' @export
simulate_max <- function(n, p, par1, par2, con=NULL) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  rdist1 <- function(n) rnorm(n, log10(par1[1]), log10(par1[2]))
  rdist2 <- function(n) rnorm(n, log10(par2[1]), log10(par2[2])) 

  rtotal <- function(n) rmix(n, p, rdist1, rdist2)

  if (is.null(con)) {
    out <- rtotal(n)
  } else {
    out <- rcon(n, rtotal, con)
  }

  return(10**out)
}

#-------------------------------------------------------------------------
#' Simulating concentration changes 
#'
#' Simulates relative (fractional) changes in metabolite concentrations using
#' a mixture of beta distributions.
#'
#' @param n Number of samples to draw.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#' @param par1 Parameters (alpha, beta) for first distribution.
#' @param par2 Parameters (alpha, beta) for second distribution.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples.
#
#' @return A vector of relative concentration values.
#'
#' @examples
#' out <- simulate_rel_change(10000, 0.7, c(2, 5), c(0.5, 0.5), c(0.1, 1))
#' hist(out, 20, probability = TRUE, 
#'      main = '', xlab = 'Fractional change in concentration')
#' @export
simulate_rel_change <- function(n, p, par1, par2, con=NULL) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  rdist1 <- function(n) rbeta(n, par1[1], par1[2])
  rdist2 <- function(n) rbeta(n, par2[1], par2[2]) 

  rtotal <- function(n) rmix(n, p, rdist1, rdist2)

  if (is.null(con)) {
    out <- rtotal(n)
  } else {
    out <- rcon(n, rtotal, con)
  }

  return(out)
}

#-------------------------------------------------------------------------
#' Simulating measurement variability 
#'
#' Simulates relative standard deviations (coefficients of variation) using
#' a mixture of 2 normal distributions (one subpopulation of less variable
#' measurements and on subpopulation of more variable ones).
#'
#' @param n Number of samples to draw.
#' @param p Proportion of samples to draw using \code{par1} vs. \code{par2}.
#' @param par1 Parameters (mean, sd) for first distribution.
#' @param par2 Parameters (mean, sd) for second distribution.
#' @param con Optional vector of lower and upper constraints on generated 
#'            samples. By default, relative standard deviations are limited
#'            to between 0 and 1.
#
#' @return A vector of relative standard deviation values.
#'
#' @examples
#' out <- simulate_rel_sd(10000, 0.7, c(.04, .02), c(0.11, 0.02), c(0, 0.20))
#' hist(out, 20, probability = TRUE, 
#'      main = '', xlab = 'Relative standard deviation')
#' @export
simulate_rel_sd <- function(n, p, par1, par2, con=c(0,1)) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  rdist1 <- function(n) rnorm(n, par1[1], par1[2])
  rdist2 <- function(n) rnorm(n, par2[1], par2[2]) 

  rtotal <- function(n) rmix(n, p, rdist1, rdist2)

  if (is.null(con)) {
    out <- rtotal(n)
  } else {
    out <- rcon(n, rtotal, con)
  }

  return(out)
}

#=========================================================================>
# Trend simulation

#-------------------------------------------------------------------------
#' Simulating decreasing trends
#'
#' Simulates decreasing trends by mixing two sigmoid curve families. Sigmoid
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function.
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max) for one set of 
#'        trends where y = (1 + exp((x-a)/b))**(-1)
#' @param par2 Parameter vector (a_min, a_max, b_min, b_max) for second set of 
#'        trends where y = (1 + exp((x-a)/b))**(-1)
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(0.2, 0.6, 0.10, 0.18)
#' par2 <- c(0.6, 0.9, 0.10, 0.18)
#' trends <- simulate_decreasing(1000, 100, 0.05, par1, par2)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_decreasing <- function(n, n.samples, p, par1, par2) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  out <- matrix(NA, nrow=n, ncol=n.samples)

  n1 <- floor(n*p)
  n2 <- n - n1

  index <- sample(1:n)
  index1 <- index[1:n1]
  index2 <- index[(n1 + 1):n]

  # Generating different sets of parameters
  par <- matrix(NA, nrow=n, ncol=2)
  par[index1, ] <- cbind(runif(n1, par1[1], par1[2]),
                         runif(n1, par1[3], par1[4]))  
  par[index2, ] <- cbind(runif(n2, par2[1], par2[2]),
                         runif(n2, par2[3], par2[4]))  

  # Defining function
  x <- seq(0, 1, length.out=n.samples)
  f_dec <- function(par) {
    y <- (1 + exp((x-par[1])/par[2]))**(-1)
    return(stretch(y, c(0, 1)))
  }

  # Generating output
  for (i in 1:n) {
    out[i, ] <- f_dec(par[i, ])
  }

  return(as.data.frame(out))
}

#-------------------------------------------------------------------------
#' Simulating increasing trends
#'
#' Simulates increasing trends by mixing two sigmoid curve families. Sigmoid
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function.
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max) for one set of 
#'        trends where y = 1 - (1 + exp((x-a)/b))**(-1)
#' @param par2 Parameter vector (a_min, a_max, b_min, b_max) for second set of 
#'        trends where y = 1 - (1 + exp((x-a)/b))**(-1)
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(0.045, 0.055, 0.2, 0.4)
#' par2 <- c(0.945, 0.955, 0.1, 0.3)
#' trends <- simulate_increasing(1000, 100, 0.15, par1, par2)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_increasing <- function(n, n.samples, p, par1, par2) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  out <- matrix(NA, nrow=n, ncol=n.samples)

  n1 <- floor(n*p)
  n2 <- n - n1

  index <- sample(1:n)
  index1 <- index[1:n1]
  index2 <- index[(n1 + 1):n]

  # Generating different sets of parameters
  par <- matrix(NA, nrow=n, ncol=2)
  par[index1, ] <- cbind(runif(n1, par1[1], par1[2]),
                         runif(n1, par1[3], par1[4]))  
  par[index2, ] <- cbind(runif(n2, par2[1], par2[2]),
                         runif(n2, par2[3], par2[4]))  

  # Defining function
  x <- seq(0, 1, length.out=n.samples)
  f_inc <- function(par) {
    y <- 1 - (1 + exp((x-par[1])/par[2]))**(-1)
    return(stretch(y, c(0, 1)))
  }

  # Generating output
  for (i in 1:n) {
    out[i, ] <- f_inc(par[i, ])
  }

  return(as.data.frame(out))
}

#-------------------------------------------------------------------------
#' Simulating concave trends
#'
#' Simulates concave trends by mixing trimmed and scaled beta distributions. 
#' parameters are sampled from uniform distributions with boundaries provided
#' as input to this function.
#'
#' @param n Number of trends to generate.
#' @param n.samples Number of timepoints within each trend.
#' @param p Proportion of trends to generate using \code{par1} vs. \code{par2}.
#' @param par1 Parameter vector (a_min, a_max, b_min, b_max, c_min, c_max, 
#'        d_min, d_max) for one set of trends where y = x**(a-1)*(1-x)**(b-1)
#'        constrained to min(x) = c, max(x) = d
#' @param par2 Parameter vector (a_min, a_max, b_min, b_max, c_min, c_max, 
#'        d_min, d_max) for one set of trends where y = x**(a-1)*(1-x)**(b-1)
#'        constrained to min(x) = c, max(x) = d
#
#' @return A dataframe with each row representing a metabolite concentration
#'         time-course scaled between 0 and 1. The corresponding x variables 
#'         are assumed to be equally spaced between 0 and 1 i.e. 
#'         x <- seq(0, 1, length.out=n.samples).    
#'
#' @examples
#' par1 <- c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9)
#' par2 <- c(3.5, 4.5, 2.5, 3.5, 0.2, 0.4, 0.7, 0.8)
#' trends <- simulate_concave(1000, 100, 0.75, par1, par2)
#'
#' # Conversion for plotting
#' y_mat <- t(as.matrix(trends))
#' x <- seq(0, 1, length.out=100)
#' matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
#'         xlab = 'Relative culturing time', ylab = 'Relative concentration')
#' @export
simulate_concave <- function(n, n.samples, p, par1, par2) {

  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  out <- matrix(NA, nrow=n, ncol=n.samples)

  n1 <- floor(n*p)
  n2 <- n - n1

  index <- sample(1:n)
  index1 <- index[1:n1]
  index2 <- index[(n1 + 1):n]

  # Generating different sets of parameters
  par <- matrix(NA, nrow=n, ncol=4)
  par[index1, ] <- cbind(runif(n1, par1[1], par1[2]),
                         runif(n1, par1[3], par1[4]), 
                         runif(n1, par1[5], par1[6]),
                         runif(n1, par1[7], par1[8]))
  par[index2, ] <- cbind(runif(n2, par2[1], par2[2]),
                         runif(n2, par2[3], par2[4]),  
                         runif(n2, par2[5], par2[6]),
                         runif(n2, par2[7], par2[8]))

  # Defining function
  x.ini <- seq(0, 1, length.out=n.samples*100)
  x.fin <- seq(0, 1, length.out=n.samples)

  f_cnc <- function(par) {
    y <- x.ini**(par[1]-1)*(1-x.ini)**(par[2]-1)

    logic <- (x.ini > par[3]) & (x.ini < par[4])
    y <- y[logic]
    x <- x.ini[logic]

    f_spline <- splinefun(stretch(x, c(0,1)), stretch(y, c(0,1)))
    out <- f_spline(x=x.fin)
    return(out)
  }

  # Generating output
  for (i in 1:n) {
    out[i, ] <- f_cnc(par[i, ])
  }

  return(as.data.frame(out))
} 

#-------------------------------------------------------------------------
#' Simulating linear trends
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

  return(as.data.frame(out))
} 



