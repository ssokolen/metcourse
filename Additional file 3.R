# Additional file 3 for BMC Systems Biology manuscript:
# Identification of internal standard addition error in 1H-NMR time-course data
# Stanislav Sokolenko, Marc G Aucoin
#
# This file contains the script used to analyze loess model performance.
# Simulation results are stored in an SQLite database.
#
# Note, this simulation took approximately 6 hours on a virtual machine
# running in a mid-range laptop.

library(dplyr)
library(ggplot2)
library(mgcv)
library(reshape2)
library(RSQLite)

# Local
source('Additional file 1.R')
source('Additional file 2.R')

#=========================================================================>
# Global options

set.seed(1111)

# Constants

n.sim <- 1000

p.d <- 0.34
p.i <- 0.30
p.c <- 0.05

bias <- 0.05

printout <- TRUE

# Changing parameters

cmp.range <- list(20, 40, 60)
obs.range <- list(10, 15, 20)
mdl.range <- list('loess-0.50', 'loess-1.0', 'loess-1.5', 'loess-2.0')

param <- expand.grid(n.cmp=cmp.range, n.obs=obs.range, mdl=mdl.range)

gen_loess <- function(span) {

  f_smooth <- function(time, concentration) {
    if (all(is.na(concentration))) return(concentration)

    d <- data.frame(time=time, concentration=concentration)
    model <- loess(concentration ~ time, data=d, span=span)
    fit <- as.vector(predict(model, d))

    return(fit)
  }

  return(f_smooth)
}

#=========================================================================>
# Simulating timecourse distributions 

for (i in 1:nrow(param)) {

  n.cmp <- param$n.cmp[[i]]
  n.obs <- param$n.obs[[i]]
  mdl <- param$mdl[[i]]
  span <- as.numeric(strsplit(mdl, '-')[[1]][2])
  f_smooth <- gen_loess(span)

  if (printout) {
    msg <- 'Generating %i simulations, n.cmpd=%i, n.obs=%i, mdl=%s\n'
    cat(sprintf(msg, n.sim, n.cmp, n.obs, mdl)) 
  }

  n.d <- ceiling(n.cmp*p.d)
  n.i <- ceiling(n.cmp*p.i)
  n.c <- ceiling(n.cmp*p.c)
  n.l <- n.cmp - n.d - n.i - n.c

  simulations <- sprintf('S%04d', 1:n.sim)
  observations <- sprintf('S%02d', 1:n.obs)

  bias <- 0.05

  #-----------------------------------------------------------------------
  # Simulating maximums

  con <- c(-2, 2)

  # Decreasing
  n <- n.sim*n.d
  p <- 0.2
  par1 <- c(-0.4, 0.35)
  par2 <- c(0.8, 0.35)
  maximums <- 10**(simulate_log10_max(n, p, par1, par2, con))

  max.d <- data.frame(maximum=maximums)
  max.d$simulation <- rep(simulations, each=n.d)
  max.d$classification <- 'dec'

  # Increasing
  n <- n.sim*n.i
  p <- 0.7
  par1 <- c(-0.8, 0.35)
  par2 <- c(0.4, 0.35)
  maximums <- 10**(simulate_log10_max(n, p, par1, par2, con))

  max.i <- data.frame(maximum=maximums)
  max.i$simulation <- rep(simulations, each=n.i)
  max.i$classification <- 'inc'

  # Concave
  n <- n.sim*n.c
  p <- 0.3
  par1 <- c(-0.8, 0.35)
  par2 <- c(0.4, 0.35)
  maximums <- 10**(simulate_log10_max(n, p, par1, par2, con))

  max.c <- data.frame(maximum=maximums)
  max.c$simulation <- rep(simulations, each=n.c)
  max.c$classification <- 'cnc'

  # Linear
  n <- n.sim*n.l
  p <- 0.3
  par1 <- c(-0.4, 0.35)
  par2 <- c(0.8, 0.35)
  maximums <- 10**(simulate_log10_max(n, p, par1, par2, con))

  max.l <- data.frame(maximum=maximums)
  max.l$simulation <- rep(simulations, each=n.l)
  max.l$classification <- 'lin'

  # Combining
  max.combined <- rbind(max.d, max.i, max.c, max.l) %>%
                    group_by(simulation) %>%
                    mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                    arrange(simulation, compound)


  #-----------------------------------------------------------------------
  # Generating percent changes

  # Decreasing
  n <- n.sim*n.d
  p <- 0.7
  par1 <- c(2, 6)
  par2 <- c(0.5, 0.5)
  con <- c(0.1, 1)
  change <- simulate_frac_change(n, p, par1, par2, con)

  change.d <- data.frame(change=change)
  change.d$simulation <- rep(simulations, each=n.d)
  change.d$classification <- 'dec'

  # Increasing
  change.i <- data.frame(change=NA, classification='inc',
                         simulation=rep(simulations, each=n.i))

  # Concave
  change.c <- data.frame(change=NA, classification='cnc',
                         simulation=rep(simulations, each=n.c))

  # Linear
  n <- n.sim*n.l
  change.l <- data.frame(change=runif(n, 0, 0.1), classification='lin',
                         simulation=rep(simulations, each=n.l))
                         
  # Combining
  change.combined <- rbind(change.d, change.i, change.c, change.l) %>%
                       group_by(simulation) %>%
                       mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                       arrange(simulation, compound)

  #-----------------------------------------------------------------------
  # Combining parameters and adding standard deviations

  parameters <- left_join(max.combined, change.combined, 
                          by=c('simulation', 'classification', 'compound'))

  n <- n.sim*n.cmp

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
  # Generating trends

  # Decreasing
  n <- n.sim*n.d
  p <- 0.05
  par1 <- c(0.20, 0.60, 0.10, 0.18)
  par2 <- c(0.60, 0.90, 0.10, 0.18)

  trend.d <- simulate_decreasing(n, n.obs, p, par1, par2)
  colnames(trend.d) <- observations

  trend.d$simulation <- rep(simulations, each=n.d)
  trend.d$classification <- 'dec'

  # Increasing
  n <- n.sim*n.i
  p <- 0.15
  par1 <- c(0.045, 0.055, 0.20, 0.40)
  par2 <- c(0.945, 0.955, 0.10, 0.30)

  trend.i <- simulate_increasing(n, n.obs, p, par1, par2)
  colnames(trend.i) <- observations

  trend.i$simulation <- rep(simulations, each=n.i)
  trend.i$classification <- 'inc'

  # Concave
  n <- n.sim*n.c
  par1 <- c(3.5, 4.5, 2.5, 3.5)
  par2 <- c(0.0, 0.2, 0.8, 0.9)

  trend.c <- simulate_concave(n, n.obs, par1, par2)
  colnames(trend.c) <- observations

  trend.c$simulation <- rep(simulations, each=n.c)
  trend.c$classification <- 'cnc'

  # Linear
  n <- n.sim*n.l
  p <- 0.6

  trend.l <- simulate_linear(n, n.obs, p)
  colnames(trend.l) <- observations

  trend.l$simulation <- rep(simulations, each=n.l)
  trend.l$classification <- 'lin'

  # Combining
  trend.combined <- rbind(trend.d, trend.i, trend.c, trend.l) %>%
                      group_by(simulation) %>%
                      mutate(compound=sprintf('C%02d', 1:n.cmp)) %>%
                      arrange(simulation, compound)

  #-----------------------------------------------------------------------
  # Combining parameters with trends, scaling trends, and adding noise

  # Joining, melting, and adding x values
  ids <- c('simulation', 'classification', 'compound')
  simulation <- left_join(parameters, trend.combined, by=ids) %>%
                melt(measure.vars=observations, 
                     variable.name='observation', value.name='trend') %>%
                arrange(simulation, compound, observation) %>%
                group_by(simulation, compound) %>%
                mutate(x=seq(0, 1, length.out=n.obs))

  # Scaling trends and adding noise
  simulation <- simulation %>%
                  group_by(simulation, compound) %>%
                  mutate(y=stretch(trend, c(minimum[1], maximum[1])),
                         deviation=rnorm(n(), 0, median(y)*st.dev[1]),
                         y=y + deviation) %>%
                  arrange(simulation, compound)

  simulation$y[simulation$y < 0] <- 0 

  #-----------------------------------------------------------------------
  # Adding bias

  # Sampling half of the simulations
  sim.biased <- sample(simulations, floor(n.sim/2))

  # Sampling one observation for each simulation
  obs.biased <- sample(observations, floor(n.sim/2), replace=TRUE)

  # Picking out appropriate rows
  combination <- paste(simulation$simulation, simulation$observation)
  logic <- combination %in% paste(sim.biased, obs.biased)

  # Adding bias
  simulation$y[logic] <- (1 + bias)*simulation$y[logic]

  # Marking biased results
  simulation$bias <- 'no'
  simulation$bias[logic] <- 'yes' 

  #-----------------------------------------------------------------------
  # Smoothing generated trends and assessing median deviation

  deviations <- simulation %>%
                  group_by(simulation, compound) %>%
                  mutate(fit=f_smooth(x, y),
                         fit.deviation=(y - fit)/y)

  # Calculating median deviation quantiles
  quants <- seq(0.05, 0.95, by=0.05)
  medians <- deviations %>%
               group_by(simulation, observation, bias) %>%
               summarize(median.deviation=
                 median(fit.deviation, na.rm=TRUE)) %>% 
               group_by(observation, bias) %>%
               do(data.frame(quantile=quants,
                             deviation=quantile(.$median.deviation, quants)))

  # Calculating corrections
  deviations <- deviations %>%
                  group_by(simulation) %>%
                  mutate(correction=correct_relative_deviation(
                           x, y, compound, max.iter=10, min.deviation=0.03)/
                           y)
  deviations$correction[is.na(deviations$correction)] <- 0

  #-----------------------------------------------------------------------
  # Formatting output to database
  drv <- dbDriver('SQLite')
  con <- dbConnect(drv, 'loess_data.db')

  # Parameters
  parameters$n.cmp <- n.cmp
  parameters$n.obs <- n.obs
  parameters$mdl <- mdl
  columns <- c('simulation', 'compound', 'n.cmp', 'n.obs', 'mdl',
               'classification', 'minimum', 'maximum', 'change', 'st.dev')  
  parameters <- parameters[, columns]

  dbWriteTable(con, 'parameters', as.data.frame(parameters), append=TRUE)

  # Trends
  deviations$n.cmp <- n.cmp
  deviations$n.obs <- n.obs
  deviations$mdl <- mdl
  columns <- c('simulation', 'compound', 'n.cmp', 'n.obs', 'mdl',
               'observation', 'x', 'y', 'fit', 'deviation', 'bias', 
               'correction')  
  deviations <- deviations[, columns]

  dbWriteTable(con, 'trends', as.data.frame(deviations), append=TRUE)

  # Medians
  medians$n.cmp <- n.cmp
  medians$n.obs <- n.obs
  medians$mdl <- mdl
  columns <- c('observation', 'bias', 'n.cmp', 'n.obs', 'mdl', 
               'deviation', 'quantile')  
  medians <- medians[, columns]

  dbWriteTable(con, 'medians', as.data.frame(medians), append=TRUE)

  dbDisconnect(con)
}
