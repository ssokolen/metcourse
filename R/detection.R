# Functions for detecting systematic bias

#========================================================================>
# Smoothing functions used for correction

#------------------------------------------------------------------------
#' GAM smoothing
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
#' LOESS smoothing
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
# Detection function

#------------------------------------------------------------------------
#' Detect systematic bias in metabolomic time-course data
#'
#' Identifies systematic deviations in metabolic data using a B-spline fit.
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
#' @param min.deviation Smallest median relative deviation to identify as a 
#'                      bias. If not supplied, it will be estimated as 50% of the underlying noise in the data.
#' @param ... Arguments to be passed into \code{f_smooth}, (such as degree and number of knots). 
#'
#' @return A dataframe corresponding to the time points with systematic error.
#'
#' @examples 
#' 
#' # Using previously simulated data 40 metabolic trends with 10 time points
#' data(timecourse)
#'
#' # Adding an error of 5% at sample 4
#' logic <- timecourse$sample == 4
#' timecourse$concentration[logic] <- timecourse$concentration[logic] * 1.05 
#'
#' # Estimating
#' error <- correct_rel_bias(timecourse$time, 
#'                                          timecourse$concentration,
#'                                          timecourse$metabolite)
#'
#'
#' @export
detect_rel_bias <- function(time, concentration, metabolite, min.deviation = NULL, ...) {

  d <- data.frame(time, concentration, 
                  metabolite, original = 1:length(time))
  
  extra_args <- list(...)

  # Passing in extra parameters into smoothing
  f_smooth_dply <- function(time, concentration) {
    args <- list(time = time, concentration = concentration)
    do.call(f_smooth, c(args, extra_args))
  }

  # Detecting systematic deviation
  # Generating fit and calculating deviations
  d <- d %>%
         group_by(metabolite) %>%
         mutate(fit = f_smooth(time, concentration),
                deviation = (concentration - fit) / concentration)

  # Identifying the timepoint with large deviation 
  deviations <- d %>%
                  group_by(time) %>%
                  summarize(med = abs(median(deviation, na.rm = TRUE))) 
  
  # Estimate threshold if not supplied by user
  if (is.null(min.deviation)) {
    mean.dev <- d %>%
      group_by(metabolite) %>%
      summarize(mean.deviation = mean(abs(deviation), na.rm = TRUE)) 
    med.dev <- median(mean.dev$mean.deviation, na.rm = TRUE)
    min.deviation <- 0.5*med.dev
  }

  deviations$error <- ifelse(deviations$med > min.deviation, TRUE, FALSE)

  out <- deviations
  
  return(out)
}

