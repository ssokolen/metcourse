# metcourse

`metcourse` version 2 provides a simple algorithm for simultaneously fitting and correcting systematic deviations in metabolomic time-courses. Specifically, it deals with the case where all metabolite concentrations are overestimated or underestimated by a constant percentage within a single sample. This can occur from inaccurate addition or processing of an internal standard as well as a result of variable solvent levels when dealing with intracellular metabolite extraction or biofluids. The correction method was validated by simulating realistic metabolic time-courses, with the simulation framework also provided in this package.

More information on this version is available in the manuscript "A comprehensive model for separating systematic bias and noise in metabolomic timecourse data -- A nonlinear B-spline mixed effect approach". An older version of the software is also available based on the manuscript "A correction method for systematic error in 1H-NMR time-course data validated through stochastic cell culture simulation" ([open access link](http://www.biomedcentral.com/1752-0509/9/51) - doi:10.1186/s12918-015-0197-4).

# Original code

Following submission of the manuscript, the original scripts were refactored into this package with some changes to enhance usability. The original code is still available under "R/original" in the form of "Additional files" mentioned throughout the manuscript. 

# Installation

Installation can be performed directly from github using the `devtools` package:

```R
library(devtools)
install_github('ssokolen/metcourse')
```

# Table of contents

1. [Detection](#detection)
2. [Correction](#correction)
3. [Simulation](#metabolite-time-course-simulation)
  1. [Trend shape](#trend-shape)
  2. [Maximum concentrations](#maximum-concentrations)
  3. [Minimum concentrations](#minimum-concentrations)
  4. [Measurement error](#measurement-error)
  5. [Putting it all together](#putting-it-all-together)

# Detection
The detection of systematic bias is performed by `detect_rel_bias`, which requires only three arguments -- time or sample number, metabolite concentration, and a column of metabolites or compounds corresponding to each concentration. 
```R
  ## Generating a set of metabolic trends to perform correction
 
  # Using previously simulated data 40 metabolic trends with 10 time points
  # (see Simulation section below and the `simulate_timecourse` example)
  data(timecourse)

  # Artificially adding an error of 5% at sample 4
  logic <- timecourse$sample == 4
  timecourse$concentration[logic] <- timecourse$concentration[logic] * 1.05
  error <- correct_rel_bias(timecourse$time, 
                                           timecourse$concentration,
                                           timecourse$metabolite)
```
# Correction

The correction of systematic bias is performed by `correct_rel_bias`, which requires only three arguments -- time or sample number, metabolite concentration, and a column of metabolites or compounds corresponding to each concentration. 

```R
  ## Generating a set of metabolic trends to perform correction
 
  # Using previously simulated data 40 metabolic trends with 10 time points
  # (see Simulation section below and the `simulate_timecourse` example)
  data(timecourse)

  # Artificially adding an error of 5% at sample 4
  logic <- timecourse$sample == 4
  timecourse$concentration[logic] <- timecourse$concentration[logic] * 1.05

  ## The correction itself

  output <- correct_rel_bias(timecourse$time,
                                           timecourse$concentration,
                                           timecourse$metabolite)
  timecourse$corrected <- output$fit

  # Plotting -- the original value of the corrected point is marked in red
  par(mfrow = c(8, 5), oma = c(5, 4, 1, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)
  new.time <- seq(min(timecourse$time), max(timecourse$time), length.out=100)

  for (metabolite in unique(timecourse$metabolite)) {

    logic <- timecourse$metabolite == metabolite
    d <- timecourse[logic, ]

    logic2 <- d$concentration != d$corrected

    plot(d$time, d$corrected, pch = 16, xlab = '', ylab = '',
    ylim = c(min(d$concentration), max(d$concentration)))

    smoothed <- met_smooth_gam(d$time, d$corrected,
    new.time = new.time, k = 4)
    lines(new.time, smoothed)

    points(d$time[logic2], d$concentration[logic2], pch = 16, col = 'red')
  }

  title(xlab = 'Time post inoculation (hours)',
        ylab = 'Concentration (mM)', outer = TRUE, line = 3)
```

# Simulation

The simulation of metabolite concentration time-courses is divided into 4 parameters -- the overall shape of the trend, maximum concentration, minimum concentration, and measurement variability (noise). Each step can be run separately or combined into a single simulation.

## Trend shape 

Metabolic time-course shapes are divided into decreasing, increasing, and concave. Decreasing and increasing trends are modelled by sigmoidal equations, y = (1 + exp((x-a)/b))^(-1) and y = 1 - (1 + exp((x-a)/b))^(-1), respectively. Concave trends are modelled by a beta distribution y = x^(a-1) * (1-x)^(b-1), truncated to allow different starting and final concentrations. 
 
```R
  # Parameters take the form of lower and upper bounds 
  # on the a and b coefficients in (1 + exp((x-a)/b))^(-1)
  par <- c(0.2, 0.6, 0.10, 0.18)

  # Simulating 1000 trends with 100 timepoints per trend
  trends <- simulate_decreasing(1000, 100, par)

  # Plotting all trends at once
  y_mat <- t(as.matrix(trends))
  x <- seq(0, 1, length.out=100)
  matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
          xlab = 'Relative culturing time', ylab = 'Relative concentration')

  # Trends can also be constructed from a mixture of two sets of parameters
  par1 <- c(0.2, 0.6, 0.10, 0.18)

  # Modifying the bounds on coefficient a
  par2 <- c(0.6, 0.9, 0.10, 0.18)

  # Only 0.05 of the trends will be generated using par1
  trends <- simulate_decreasing(1000, 100, par1, par2, 0.05)

  # Plotting
  y_mat <- t(as.matrix(trends))
  x <- seq(0, 1, length.out=100)
  matplot(x, y_mat, type = 'l', lty = 1, lwd = 4, col = grey(0, 0.05),
          xlab = 'Relative culturing time', ylab = 'Relative concentration')
```

See examples for `simulate_increasing` and `simulate_concave` i.e. `example(simulate_increasing)` and `example(simulate_concave)`.

## Maximum concentrations

Maximum metabolite concentrations are generated from a normal distribution of concentration logarithms.

```R
  # Parameters take the form of mean and standard deviation values
  # in units of concentration (base-10 logarithms are taken internally)
  par <- c(7, 2)
  out <- simulate_max(10000, par)

  # Formatting and plotting output on logarithmic scale
  out <- log10(out)
  labels <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50)
  hist(out, 20, axes = FALSE, probability = TRUE, 
       main = '', xlab = 'Maximum metabolite concentration')
  axis(side = 2)
  axis(at = log10(labels), labels = labels, side = 1)

  # A mixture of two distributions can be used to account for a large
  # number of low concentration metabolites. 
  par1 <- c(7, 2)
  par2 <- c(0.5, 2)

  # Only 0.3 of the maximum concentration are sampled from the
  # higher concentration distribution.
  out <- simulate_max(10000, par1, par2, 0.3)

  # Formatting and plotting output on logarithmic scale
  out <- log10(out)
  labels <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50)
  hist(out, 20, axes = FALSE, probability = TRUE, 
       main = '', xlab = 'Maximum metabolite concentration')
  axis(side = 2)
  axis(at = log10(labels), labels = labels, side = 1)

  # Constraints can be set on upper and lower values
  # (while maintaining the number of maximum concentrations generated)
  out <- simulate_max(10000, par1, par2, 0.3, con=c(0.1, 20))

  # Formatting and plotting output on logarithmic scale
  out <- log10(out)
  labels <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50)
  hist(out, 20, axes = FALSE, probability = TRUE, 
       main = '', xlab = 'Maximum metabolite concentration')
  axis(side = 2)
  axis(at = log10(labels), labels = labels, side = 1)
```

## Minimum concentrations

To avoid correlations between minimum and maximum concentrations, minimum concentration are calculated by generating a relative (fractional) change from maximum concentration. A mixture of two beta distributions gives a fine degree of control over the distribution.

```R
  # The syntax for simulate_rel_change is very similar to simulate_max
  
  # First set of parameters to favour changes in the 25% range
  par1 <- c(2, 5)

  # Second set of parameters to add high probability of either no
  # change, or very large change (optional)
  par2 <- c(0.5, 0.5)

  # (Constraints can also be used to set minimum/maximum changes)

  # Simulating
  out <- simulate_rel_change(10000, par1, par2, 0.7)

  # Plotting
  hist(out, 20, probability = TRUE,
       main = 'Relative change', xlab = 'Fractional change in concentration')

  # Repeating with only one set of parameters and with constraints
  out <- simulate_rel_change(10000, par1, con = c(0.1, 0.5))

  # Plotting
  hist(out, 20, probability = TRUE,
       main = 'Relative change', xlab = 'Fractional change in concentration')
```

## Measurement error

The simulation of measurement error follows the same pattern as above, with the distribution of relative standard deviations modelled as a mixture of normal distributions.

```R
  # First set of parameters for high precision measurements (sd ~ 5%)
  par1 <- c(0.05, 0.02)

  # Second set of parameters for lower precision measurements (sd ~ 10%)
  par2 <- c(0.10, 0.04)

  # Constraints are important to prevent negative or impossibly low variability
  con <- c(0.01, 1)

  # Simulating
  out <- simulate_rel_sd(10000, par1, par2, 0.7, con)

  # Plotting
  hist(out, 20, probability = TRUE,
       main = 'Measurement variability', xlab = 'Fractional standard deviation')
```

## Putting it all together

The above functions can be combined manually to generate a set of trends. Alternatively, `simulate_timecourse` simplifies the process, but a large number of parameters are still required. 

```R
  # This example makes use of tidyr and dplyr packages
  library(dplyr)
  library(tidyr)

  # Making up metabolite names
  metabolites <- paste('met', 1:6, sep = '')

  # Generating a set of 6 decreasing trends with 20 timepoints
  trends <- simulate_decreasing(n = 6, n.samples = 20, 
                                par1 = c(0.2, 0.6, 0.10, 0.18))

  # Changing column names of samples
  colnames(trends) <- paste('sample', 1:20, sep = '')

  # Adding metabolite column
  trends$metabolite <- metabolites

  # Initializing data frame of parameters
  parameters <- data.frame(metabolite = metabolites, stringsAsFactors = FALSE)

  # Adding maximums
  parameters$maximum <- simulate_max(n = 6, par1 = c(7, 2), con = c(0, 50))

  # Adding minimums 
  parameters$change <- simulate_rel_change(n = 6, par1 = c(5, 2))
  parameters$minimum <- (1 - parameters$change) * parameters$maximum

  # Adding measurement error
  parameters$sdv <- simulate_rel_sd(n = 6, par1 = c(0.02, 0.02), 
                                    con = c(0.01, 1))

  # Combining
  combined <- left_join(trends, parameters, by = 'metabolite')  

  # Converting into "long" format
  combined <- gather(combined, sample, conc, sample1:sample20)

  # Scaling trends
  combined <- combined %>%
                group_by(metabolite) %>%
                mutate(conc = conc * (maximum - minimum) + minimum,
                       error = rnorm(20, 0, mean(conc) * sdv[1]),
                       conc = conc + error,
                       sample = as.numeric(gsub('sample', '', sample)))

  # Plotting
  par(mfrow = c(3, 2), oma = c(5, 4, 1, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)

  for (metabolite in unique(combined$metabolite)) {
    logic <- combined$metabolite == metabolite
    plot(combined$sample[logic], combined$conc[logic],
    xlab='', ylab='')
  }

  title(xlab = 'Sample', ylab = 'Concentration', outer = TRUE, line = 3)
```

The following example simulates a full metabolic profile from a suspension cell culture (modelled on the growth of insect cells). The number of metabolites is limited to 10 for easier plotting.

```R
  # All the parameters are defined at once in a single list. Global parameters
  # in the form of `par1.max` are used in the generation of maximum values for
  # all trend types. These can be overriden for specific types of trends by
  # appending the trend type to the parameter e.g. `par1.max.increasing`.

  param <- list(
    # Maximum concentrations are the same for every trend type
    p.max = 0.3,
    par1.max = c(7, 2),
    par2.max = c(0.5, 2),
    con.max = c(0, 50),

    # Global change parameters are near 100% for increasing/concave trends
    par1.change = c(5, 0.1),
    con.change = c(0.5, 1),

    # Decreasing trends can have a wide variety of changes
    p.change.decreasing = 0.7,
    par1.change.decreasing = c(2, 5),
    par2.change.decreasing = c(0.5, 0.5),
    con.change.decreasing = c(0.1, 1),

    # Linear trends are characterized by relatively small changes
    par1.change.linear = c(1, 5),
    con.change.linear = c(0, 0.1),

    # Measurement error is the same for every trend type (but no more than 20%)
    p.sd = 0.7,
    par1.sd = c(0.04, 0.02),
    par2.sd = c(0.11, 0.02),
    con.sd = c(0, 0.20),

    # Decreasing trend specification
    p.trend.decreasing = 0.05,
    par1.trend.decreasing = c(0.2, 0.6, 0.10, 0.18),
    par2.trend.decreasing = c(0.6, 0.9, 0.10, 0.18),

    # Increasing trend specification
    p.trend.increasing = 0.15,
    par1.trend.increasing = c(0.045, 0.055, 0.2, 0.4),
    par2.trend.increasing = c(0.945, 0.955, 0.1, 0.3),

    # Concave trend specification
    par1.trend.concave = c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9),

    # Linear trends are equaly split between increasing and decreasing
    p.trend.linear = 0.5
  )

  # Generating trends
  timecourse <- simulate_timecourse(10, c(0.3, 0.3, 0.3), param)

  # Plotting
  par(mfrow = c(5, 2), oma = c(5, 4, 1, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)

  for (metabolite in unique(timecourse$metabolite)) {
    logic <- timecourse$metabolite == metabolite
    plot(timecourse$sample[logic], timecourse$concentration[logic],
    xlab='', ylab='')
  }

  title(xlab = 'Sample', ylab = 'Concentration', outer = TRUE, line = 3)
```


