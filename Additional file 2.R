# Additional file 2 for BMC Systems Biology manuscript:
# Identification of internal standard addition error in 1H-NMR time-course data
# Stanislav Sokolenko, Marc G Aucoin
#
# This file contains functions used to simulate realistic cell culture
# time course trends.
#
# An example script can be found at the end of the file.

#=========================================================================>
# Helper functions

#-------------------------------------------------------------------------
rcon <- function(n, rdist, con, max.iter=100) {
  # Resamples supplied distribution to ensure generated values fall
  # within given constraints.
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of samples to draw
  # rdist: function(n)
  #   User-defined function for generating random samples
  # con: numeric
  #   Lower and upper constraints on generated samples
  # max.iter: integer
  #   Maximum number of iterations before aborting

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

#-------------------------------------------------------------------------
rmix <- function(n, p, rdist1, rdist2) {
  # Generates a random sample with samples coming from rdist1 p fraction
  # of the time and rdist2 1-p fraction of the time.
  #
  # Parameters
  # ----------
  # n: integer
  #   Number of samples to draw
  # p: numeric
  #   Approximate proportion of samples to draw from rdist1
  # rdist1: function(n)
  #   User-defined function for generating random samples
  # rdist2: function(n)
  #   User-defined function for generating random samples

  split <- runif(n)
  n1 <- sum(split < p)
  n2 <- n - n1

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

#=========================================================================>
# An example script to generate the trends from a single experiment 

if (FALSE) {

  # Setting seed for random generation
  set.seed(1113)

  # Packages required to run example
  check.dplyr <- require(dplyr)
  check.ggplot <- require(ggplot2)
  check.reshape <- require(reshape2)

  if (!all(check.dplyr, check.ggplot, check.reshape)) {
    msg <- 'The packages dplyr, ggplot2, and reshape2 are required for example'
    stop(msg)
  }

  # Setting parameters 
  n.cmp <- 40
  n.obs <- 10
  n.trends <- n.cmp * n.obs

  # Trend percentages
  p.decreasing <- 0.34
  p.increasing <- 0.30
  p.concave <- 0.05

  # Converting to numbers
  n.decreasing <- ceiling(n.cmp*p.decreasing)
  n.increasing <- ceiling(n.cmp*p.increasing)
  n.concave <- ceiling(n.cmp*p.concave)
  n.linear <- n.cmp - n.decreasing - n.increasing - n.concave

  #-----------------------------------------------------------------------
  # Generating trends

  # Defining sample column names
  observations <- sprintf('S%02d', 1:n.obs)

  # Decreasing
  p <- 0.05
  par1 <- c(0.20, 0.60, 0.10, 0.18)
  par2 <- c(0.60, 0.90, 0.10, 0.18)

  trend.decreasing <- simulate_decreasing(n.decreasing, n.obs, p, par1, par2)
  colnames(trend.decreasing) <- observations
  trend.decreasing$classification <- 'decreasing'

  # Increasing
  p <- 0.15
  par1 <- c(0.045, 0.055, 0.20, 0.40)
  par2 <- c(0.945, 0.955, 0.10, 0.30)

  trend.increasing <- simulate_increasing(n.increasing, n.obs, p, par1, par2)
  colnames(trend.increasing) <- observations
  trend.increasing$classification <- 'increasing'

  # Concave
  par1 <- c(3.5, 4.5, 2.5, 3.5)
  par2 <- c(0.0, 0.2, 0.8, 0.9)

  trend.concave <- simulate_concave(n.concave, n.obs, par1, par2)
  colnames(trend.concave) <- observations
  trend.concave$classification <- 'concave'

  # Linear
  p <- 0.6

  trend.linear <- simulate_linear(n.linear, n.obs, p)
  colnames(trend.linear) <- observations
  trend.linear$classification <- 'linear'

  # Combining all trends together
  trend.combined <- rbind(trend.decreasing, trend.increasing, 
                          trend.concave, trend.linear) %>%
                      mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                      arrange(compound)

  #-----------------------------------------------------------------------
  # Simulating maximums

  # Setting log10(concentration) constraints
  con <- c(-2, 2)

  # Setting bimodal parameters for all trends
  par1 <- c(-0.4, 0.35)
  par2 <- c(0.8, 0.35)

  # Decreasing
  p <- 0.2
  maximums <- 10**(simulate_log10_max(n.decreasing, p, par1, par2, con))

  max.decreasing <- data.frame(maximum=maximums)
  max.decreasing$classification <- 'decreasing'

  # Increasing
  p <- 0.7
  maximums <- 10**(simulate_log10_max(n.increasing, p, par1, par2, con))

  max.increasing <- data.frame(maximum=maximums)
  max.increasing$classification <- 'increasing'

  # Concave
  p <- 0.3
  maximums <- 10**(simulate_log10_max(n.concave, p, par1, par2, con))

  max.concave <- data.frame(maximum=maximums)
  max.concave$classification <- 'concave'

  # Linear
  p <- 0.3
  maximums <- 10**(simulate_log10_max(n.linear, p, par1, par2, con))

  max.linear <- data.frame(maximum=maximums)
  max.linear$classification <- 'linear'

  # Combining all the maximums
  max.combined <- rbind(max.decreasing, max.increasing, 
                        max.concave, max.linear) %>%
                    mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                    arrange(compound)

  #-----------------------------------------------------------------------
  # Generating percent changes

  # Decreasing
  p <- 0.7
  par1 <- c(2, 6)
  par2 <- c(0.5, 0.5)
  con <- c(0.1, 1)
  change <- simulate_frac_change(n.decreasing, p, par1, par2, con)

  change.decreasing <- data.frame(change=change)
  change.decreasing$classification <- 'decreasing'

  # Increasing
  change.increasing <- data.frame(change=rep(NA, n.increasing), 
                                  classification='increasing')

  # Concave
  change.concave <- data.frame(change=rep(NA, n.concave), 
                               classification='concave')
  
  # Linear
  change.linear <- data.frame(change=runif(n.linear, 0, 0.1), 
                              classification='linear')
                         
  # Combining
  change.combined <- rbind(change.decreasing, change.increasing, 
                           change.concave, change.linear) %>%
                       mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                       arrange(compound)

  #-----------------------------------------------------------------------
  # Combining parameters and adding standard deviations

  parameters <- left_join(max.combined, change.combined, 
                          by=c('classification', 'compound'))

  p <- 0.65
  par1 <- c(.04, .02)
  par2 <- c(.11, .02)
  con <- c(.01, .20)

  parameters$st.dev <- simulate_rel_sd(n.cmp, p, par1, par2, con)

  # Calculating minimum values
  parameters <- within(parameters, {
    minimum <- maximum*(1 - change)
    minimum[is.na(change)] <- 0
  })

  #-----------------------------------------------------------------------
  # Combining parameters with trends, scaling trends, and adding noise

  # Joining, melting, and adding x values
  ids <- c('classification', 'compound')
  simulation <- left_join(parameters, trend.combined, by=ids) %>%
                melt(measure.vars=observations, 
                     variable.name='observation', value.name='trend') %>%
                arrange(compound, observation) %>%
                group_by(compound) %>%
                mutate(x=seq(0, 1, length.out=n.obs))

  # Scaling trends and adding noise
  simulation <- simulation %>%
                  group_by(compound) %>%
                  mutate(y=stretch(trend, c(minimum[1], maximum[1])),
                         deviation=rnorm(n(), 0, median(y)*st.dev[1]),
                         y=y + deviation) %>%
                  arrange(compound)

  simulation$y[simulation$y < 0] <- 0 

  #-----------------------------------------------------------------------
  # Renaming compounds in order of concentration
  compounds <- sort(unique(simulation$compound))
  new.order <- simulation %>% 
                 group_by(compound) %>%
                 summarize(maximum=maximum[1]) %>%
                 arrange(desc(maximum))
  new.order <- new.order$compound
  names(compounds) <- new.order

  simulation$compound <- compounds[simulation$compound]

  #-----------------------------------------------------------------------
  # Generating output

  p <- ggplot(simulation, aes(x=x, y=y))
  p <- p + geom_point(size=3)

  # loess smoothing used instead of gam to avoid requiring extra mgcv package
  p <- p + stat_smooth(method='loess', span=0.75, colour='black')

  p <- p + xlab('Relative time')
  p <- p + ylab('Concentration (mM)')

  p <- p + facet_wrap(~ compound, ncol=3, scale='free_y')
  p <- p + theme_bw(18)

  ggsave('example_trends.pdf', width=10, height=20, units='in')

}




    

