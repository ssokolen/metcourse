# metcourse

Scripts for simulating metabolic time-course data from suspension cell cultures and correcting systematic deviations based on the manuscript "A correction method for systematic error in 1H-NMR time-course data validated through stochastic cell culture simulation".

# Original code

Following submission of the manuscript, the original scripts were refactored into this package with some changes to enhance usability. The original code is still available under "R/original" in the form of "Additional files" mentioned throughout the manuscript. 

# Installation

Installation can be performed directly from github using the `devtools` package:

```R
library(devtools)
install_github('ssokolen/metcourse')
```

# Metabolite time-course simulation

The simulation of metabolite concentration time-courses is divided into 4 parameters -- the overall shape of the trend, maximum concentration, minimum concentration, and measurement variability (noise). Each step can be run separately or all at once.

## Trend shape 

Metabolic time-course shapes are divided into decreasing, increasing, and concave. Decreasing and increasing trends are modelled by sigmoidal equations, y = (1 + exp((x-a)/b))^(-1) and y = 1 - (1 + exp((x-a)/b))^(-1), respectively. Concave trends are modelled by a beta distribution y = x^(a-1) * (1-x)^(b-1), truncated to allow different starting and final concentrations. 
 
```R
  # Parameters take the form of lower and upper bounds 
  # on the a and b coefficients in (1 + exp((x-a)/b))^(-1)
  par <- c(0.2, 0.6, 0.10, 0.18)

  # Simulating 1000 trends with 100 timepoints per trend using
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

See examples for `simulate_increasing` and `simulate_concave`

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

