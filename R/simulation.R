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
#' # Using default distributions parameters
#' out <- rcon(10000, rnorm, c(-1, 1))
#' print(summary(out))
#' hist(out, 20)
#'
#' # Modifying mean and standard deviation
#' out <- rcon(10000, function(n) {rnorm(n, 0.5, 0.5)}, c(-1, 1))
#' print(summary(out))
#' hist(out, 20)
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
#' @param p Approximate proportion of samples to draw from rdist1 vs rdist2 
#'          (true proportion taken as the fraction of random uniform values
#'           smaller than p).
#' @param rdist1 User-defined function for generating random samples 
#'              (must take \code{n} as only argument).
#' @param rdist2 User-defined function for generating random samples 
#'              (must take \code{n} as only argument).
#' @param exact TRUE if \code{p} should be interpreted as the exact proportion.
#
#' @return A vector of \code{n} values sampled from \code{rdist1} and
#'         \code{rdis2}.
#'
#' @examples
#' rdist1 <- function(n) {rnorm(n, -1, 0.5)}
#' rdist2 <- function(n) {rnorm(n,  1, 0.5)}
#' out <- rmix(10000, 0.3, rdist1, rdist2)
#' print(summary(out))
#' hist(out, 20)
#' @export
rmix <- function(n, p, rdist1, rdist2, exact=FALSE) {
  
  if ((p < 0) | (p > 1)) {
   msg <- 'p must be between 0 and 1'
   stop(msg)
  } 

  if (exact) {
    n1 <- floor(n*p)
    n2 <- n - n1
  } else {
    split <- runif(n)
    n1 <- sum(split < p)
    n2 <- n - n1
  }

  out <- c(rdist1(n1), rdist2(n2))
  
  return(sample(out, n))
}

#-------------------------------------------------------------------------
stretch <- function(x, con) {
  # Maps x values between specified constraints.
  #
  # Parameters
  # ----------
  # x: numeric
  #   Data
  # con: numeric
  #   Vector in the form of c(min, max) specifying new constraints

  x <- (x - min(x))/(max(x) - min(x)) * (con[2] - con[1]) + con[1]
  return(x)
}

#=========================================================================>
# Parameter simulation

#-------------------------------------------------------------------------
simulate_log10_max <- function(n, p, par1, par2, con=NULL) {
  # Simulates logarithms of maximum compound concentrations using a 
  # mixture of 2 normal distributions
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of samples to generate
  # p: numeric
  #   Fraction of sample to generate using par1 vs par2
  # par1: numeric
  #   Parameter vector for first normal distribution (mean, sd)
  # par2: numeric
  #   Parameter vector for second normal distribution (mean, sd)

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

#-------------------------------------------------------------------------
simulate_frac_change <- function(n, p, par1, par2, con=NULL) {
  # Simulates fractional reduction in concentration of decreasing 
  # compounds using a mixture of 2 beta distributions
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of samples to generate
  # p: numeric
  #   Fraction of sample to generate using par1 vs par2
  # par1: numeric
  #   Parameter vector for first beta distribution (alpha, beta)
  # par2: numeric
  #   Parameter vector for second beta distribution (alpha, beta)

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
simulate_rel_sd <- function(n, p, par1, par2, con=NULL) {
  # Simulates relative standard deviation (as a fraction of mean 
  # concentration) using a mixture of 2 normal distributions.
  # [This function is currently identical to simulate_log10_max]
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of samples to generate
  # p: numeric
  #   Fraction of sample to generate using par1 vs par2
  # par1: numeric
  #   Parameter vector for first normal distribution (mean, sd)
  # par2: numeric
  #   Parameter vector for second normal distribution (mean, sd)

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
simulate_decreasing <- function(n, n.samples, p, par1, par2) {
  # Simulates decreasing trends by mixing two sigmoid curve families.
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of trends to generate 
  # n.samples: integer
  #   Number of timepoints within each trend
  # p: numeric
  #   Fraction of sample to generate using par1 vs par2
  # par1: numeric
  #   Parameter vector for one set of trends (a_min, a_max, b_min, b_max)
  #   where y = (1 + exp((x-a)/b))**(-1)
  # par2: numeric
  #   Parameter vector for second set of trends (a_min, a_max, b_min, b_max)
  #   where y = (1 + exp((x-a)/b))**(-1)

  out <- matrix(NA, nrow=n, ncol=n.samples)

  split <- runif(n)
  
  logic1 <- split < p
  n1 <- sum(logic1)

  logic2 <- split >= p
  n2 <- sum(logic2)

  # Generating different sets of parameters
  par <- matrix(NA, nrow=n, ncol=2)
  par[logic1, ] <- cbind(runif(n1, par1[1], par1[2]),
                         runif(n1, par1[3], par1[4]))  
  par[logic2, ] <- cbind(runif(n2, par2[1], par2[2]),
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
simulate_increasing <- function(n, n.samples, p, par1, par2) {
  # Simulates increasing trends by mixing two sigmoid curve families.
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of trends to generate 
  # n.samples: integer
  #   Number of timepoints within each trend
  # p: numeric
  #   Fraction of sample to generate using par1 vs par2
  # par1: numeric
  #   Parameter vector for one set of trends (a_min, a_max, b_min, b_max)
  #   where y = 1 - (1 + exp((x-a)/b))**(-1)
  # par2: numeric
  #   Parameter vector for second set of trends (a_min, a_max, b_min, b_max)
  #   where y = 1 - (1 + exp((x-a)/b))**(-1)

  out <- matrix(NA, nrow=n, ncol=n.samples)

  split <- runif(n)
  
  logic1 <- split < p
  n1 <- sum(logic1)

  logic2 <- split >= p
  n2 <- sum(logic2)

  # Generating different sets of parameters
  par <- matrix(NA, nrow=n, ncol=2)
  par[logic1, ] <- cbind(runif(n1, par1[1], par1[2]),
                         runif(n1, par1[3], par1[4]))  
  par[logic2, ] <- cbind(runif(n2, par2[1], par2[2]),
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
simulate_concave <- function(n, n.samples, par1, par2) {
  # Simulates increasing trends by trimming and scaling a beta distribution.
  # The first set of parameters determine the beta parameters while the
  # second set constrains trimming. Following trimming, the x values are
  # rescaled from 0 to 1
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of trends to generate 
  # n.samples: integer
  #   Number of timepoints within each trend
  # par1: numeric
  #   Parameter vector for beta distribution (a_min, a_max, b_min, b_max)
  #   where y = x**(a-1)*(1-x)**(b-1)
  # par2: numeric
  #   Parameter vector for truncation (c_min, c_max, d_min, d_max)
  #   where min(x) = c, max(x) = d


  out <- matrix(NA, nrow=n, ncol=n.samples)

  # Generating parameters
  par <- cbind(runif(n, par1[1], par1[2]),
               runif(n, par1[3], par1[4]),
               runif(n, par2[1], par2[2]),
               runif(n, par2[3], par2[4]))

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
simulate_linear <- function(n, n.samples, p) {
  # Simulates linear trends with p fraction decreasing and (1-p) fraction
  # increasing.
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of trends to generate 
  # n.samples: integer
  #   Number of timepoints within each trend
  # p: numeric
  #   Fraction of trends decreasing

  n1 <- ceiling(p*n)
  n2 <- n - n1

  # Defining function
  x <- seq(0, 1, length.out=n.samples)
  
  y1 <- seq(1, 0, length.out=n.samples)
  y2 <- x

  out1 <- matrix(rep(y1, each=n1), nrow=n1, ncol=n.samples)
  out2 <- matrix(rep(y2, each=n2), nrow=n2, ncol=n.samples)

  out <- rbind(out1, out2)

  return(as.data.frame(out[sample(1:n), ]))
} 

