# Functions for correcting systematic bias

#========================================================================>
# Smoothing functions used for correction

#------------------------------------------------------------------------
#' GAM Smoothing
#'
#' Wraps GAM smoothing function for use in bias correction algorithm.
#'
#' @param time Time or sample number for metabolite time-course.
#' @param concentration Metabolite concentration.
#' @param new.time Optional vector for new independent variables.
#' @param warn The \code{gam} function uses internal warning suppression that
#'             doesn't appear to work when used in conjunction with dplyr.
#'             All warning messages are suppressed by default here. Set to
#'             TRUE to unsuppress.
#' @param ... Arguments to be passed into the s() smoothing function, such as
#'            k -- the dimension of the basis used to represent the smooth term.
#'            If no additional arguments are provided, a reasonable default
#'            for k (5) is assumed.
#
#' @return A vector of smoothed concentrations.
#'
#' @examples
#' # Simulating concave curve with 10 points
#' par <- c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9)
#' concentration <- as.numeric(simulate_concave(1, 10, par))
#' 
#' # Adding noise and changing scale
#' concentration <- concentration * 18 + 2
#' abs.sd <- mean(concentration) * 0.1
#' concentration <- concentration + rnorm(10, 0, abs.sd)
#' concentration[concentration < 0] <- 0 
#' 
#' # Original sampling time
#' time <- seq(0, 9, length.out=10) * 24
#'
#' # New time for smoothing
#' new.time <- seq(0, 9, length.out=100) * 24
#' 
#' # Smoothing
#' corrected <- met_smooth_gam(time, concentration, new.time = new.time, k = 5)
#' 
#' # Plotting
#' plot(time, concentration, type='p', 
#'      xlab='Time post inoculation (hours)', ylab='Concentration (mM)')
#' lines(new.time, corrected, type='l')
#' @export
met_smooth_gam <- function(time, concentration, new.time = NULL,
                           warn = FALSE, ...) {
  
  d <- data.frame(time = time, concentration = concentration)

  if (warn) {
    f_gam <- function(...) gam(...)
  } else {
    f_gam <- function(...) suppressWarnings(gam(...))
  }

  if (length(list(...)) == 0) {
    model <- f_gam(concentration ~ s(time, bs='cr', k = 5), data = d)
  } else {
    model <- f_gam(concentration ~ s(time, bs='cr', ...), data = d)
  }

  if (! is.null(new.time)) {
    d <- data.frame(time = new.time)
  }

  fit <- as.vector(predict(model, d))

  return(fit)
}

#------------------------------------------------------------------------
#' LOESS Smoothing
#'
#' Wraps LOESS smoothing function for use in bias correction algorithm. 
#' WARNING: LOESS is not a particularly good choice for bias correction. 
#' This function is included for completeness only.
#'
#' @param time Time or sample number for metabolite time-course.
#' @param concentration Metabolite concentration.
#' @param new.time Optional vector for new independent variables.
#' @param ... Arguments to be passed into the loess() smoothing function, 
#'            such as span -- "the parameter alpha which controls the degree of
#'            smoothing". If no additional arguments are provided, a reasonable 
#'            default for span (0.75) is assumed.
#
#' @return A vector of smoothed concentrations.
#'
#' @examples
#' # Simulating concave curve with 10 points
#' par <- c(3.5, 4.5, 2.5, 3.5, 0.0, 0.2, 0.8, 0.9)
#' concentration <- as.numeric(simulate_concave(1, 10, par))
#' 
#' # Adding noise and changing scale
#' concentration <- concentration * 18 + 2
#' abs.sd <- mean(concentration) * 0.1
#' concentration <- concentration + rnorm(10, 0, abs.sd)
#' concentration[concentration < 0] <- 0 
#' 
#' # Original sample time
#' time <- seq(0, 9, length.out=10) * 24
#'
#' # New time for smoothing
#' new.time <- seq(0, 9, length.out=100) * 24
#' 
#' # Smoothing
#' corrected <- met_smooth_loess(time, concentration, 
#'                               new.time = new.time, span = 0.75)
#' 
#' # Plotting
#' plot(time, concentration, type='p', 
#'      xlab='Time post inoculation (hours)', ylab='Concentration (mM)')
#' lines(new.time, corrected, type='l')
#' @export
met_smooth_loess <- function(time, concentration, new.time = NULL, ...) {
  
  d <- data.frame(time = time, concentration = concentration)

  if (length(list(...)) == 0) {
    model <- loess(concentration ~ time, data = d, span = 0.75)
  } else {
    model <- loess(concentration ~ time, data = d, ...)
  }

  if (! is.null(new.time)) {
    d <- data.frame(time = new.time)
  }

  fit <- as.vector(predict(model, d))

  return(fit)
}

#========================================================================>
# Correction function

#------------------------------------------------------------------------
#' Correcting Systematic Bias in Metabolomic Timecourse Data
#'
#' Identifies systematic deviations in metabolic data using a smoothing fit.
#' Timepoints that have a median relative deviation above a threshold value
#' are assumed to be influenced by a measurement or methodological bias as
#' compared to the overall trends metabolite concentrations.
#'
#' @param time A vector of times or sample numbers for metabolite time-courses.
#'             Note, there should be a time value for every \code{concentration} 
#'             for every \code{metabolite} e.g. time = c(1, 2, 3, 1, 2, 3), 
#'             concentration = c(20, 15, 10, 3, 6, 9), 
#'             metabolite = c('glc', 'glc', 'glc', 'lac', 'lac', 'lac'). 
#' @param concentration A vector of metabolite concentrations.
#' @param metabolite A vector of metabolite names that correspond to \code{time}
#'                   and \code{concentration}.
#' @param f_smooth Smoothing function (set to met_smooth_gam by default).
#' @param max.iter The algorithm correct the single most deviating point at one
#'                 time (as it may influence the identification of other
#'                 deviations). \code{max.iter} controls the maximum number
#'                 of correction passes.
#' @param min.deviation Smallest median relative deviation to identify as a 
#'                      bias. 0.02 has been found to be a good default for
#'                      a diverse set of metabolic time-courses from
#'                      higher eukaryote cell culture.
#' @param ... Arguments to be passed into \code{f_smooth}, (such as k or
#'            span parameters for met_smooth_gam and met_smooth_loess
#'            respectively). 
#
#' @return A vector of corrected concentrations.
#'
#' @examples
#' set.seed(1111)
#' 
#' # Setting parameters to generate 40 metabolic trends with 10 time points
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
#' timecourse <- simulate_timecourse(40, c(0.4, 0.2, 0.2), param)
#'
#' # Adding times
#' timecourse$time <- (timecourse$sample - 1)*24
#'
#' # Adding an error of 10% at sample 4
#' logic <- timecourse$sample == 4
#' timecourse$concentration[logic] <- timecourse$concentration[logic] * 1.05 
#'
#' # Correcting
#' timecourse$corrected <- correct_rel_bias(timecourse$time, 
#'                                          timecourse$concentration,
#'                                          timecourse$metabolite)
#'
#'
#' # Plotting -- the original value of the corrected point is marked in red
#' par(mfrow = c(8, 5), oma = c(5, 4, 1, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)
#' new.time <- seq(min(timecourse$time), max(timecourse$time), length.out=100)
#' 
#' for (metabolite in unique(timecourse$metabolite)) {
#'
#'   logic <- timecourse$metabolite == metabolite
#'   d <- timecourse[logic, ]
#'   
#'   logic2 <- d$concentration != d$corrected
#'
#'   plot(d$time, d$corrected, pch = 16, xlab = '', ylab = '',  
#'        ylim = c(min(d$concentration), max(d$concentration)))
#'
#'   smoothed <- met_smooth_gam(d$time, d$corrected,
#'                              new.time = new.time, k = 5)
#'   lines(new.time, smoothed)
#'
#'   points(d$time[logic2], d$concentration[logic2], pch = 16, col = 'red')
#'
#' }
#'
#' title(xlab = 'Time post inoculation (hours)', 
#'       ylab = 'Concentration (mM)', outer = TRUE, line = 3)
#' @export
correct_rel_bias <- function(time, concentration, metabolite,
                             f_smooth = met_smooth_gam, 
                             max.iter = 20, min.deviation = 0.02, ...) {

  # Generating data frame, with a new "corrected" column
  d <- data.frame(time, concentration, corrected = concentration, 
                  metabolite, original = 1:length(time))
  
  extra_args <- list(...)

  # Passing in extra parameters into smoothing
  f_smooth_dply <- function(time, corrected) {
    args <- list(time = time, corrected = corrected)
    do.call(f_smooth, c(args, extra_args))
  }
  
  # Iteratively correcting systematic deviation
  for (i in 1:max.iter) {

    # Generating fit and calculating deviations
    d <- d %>%
           group_by(metabolite) %>%
           mutate(fit = f_smooth(time, corrected),
                  deviation = (corrected - fit) / corrected)

    # Identifying the timepoint with largest deviation 
    deviations <- d %>%
                    group_by(time) %>%
                    summarize(med = abs(median(deviation, na.rm = TRUE))) %>%
                    arrange(desc(med))

    if (deviations$med[1] < min.deviation) break

    max_dev_time <- deviations$time[1]
   
    logic <- d$time == max_dev_time
    d[logic, ] <- within(d[logic, ], {
      corrected <- corrected - median(deviation) * corrected
    })
  }

  # Correcting any changes in order
  out <- d$corrected[order(d$original)]
  
  return(out)
}

