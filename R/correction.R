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
#' @param ... Arguments to be passed into the s() smoothing function, such as
#'            k -- the dimension of the basis used to represent the smooth term.
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
met_smooth_gam <- function(time, concentration, new.time = NULL, ...) {
  
  d <- data.frame(time = time, concentration = concentration)
  model <- gam(concentration ~ s(time, bs='cr', ...), data=d)

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
#'            smoothing".
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
  
  d <- data.frame(time=time, concentration=concentration)
  model <- loess(concentration ~ time, data=d, ...)

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
#'            respectively). Note, the function attempts to identify the use of 
#'            met_smooth_gam and met_smooth_loess in order to provide reasonable
#'            default smoothing parameters based on observations of metabolic 
#'            time-courses from higher eukaryote cell culture. It is highly
#'            recommened to provide explicit smoothing parameters.
#
#' @return A vector of corrected concentrations.
#'
#' @examples
#' # TODO
#' @export
correct_rel_deviation <- function(time, concentration, metabolite,
                                  f_smooth = met_smooth_gam, 
                                  max.iter=20, min.deviation=0.02, ...) {

  # Checking f_smooth to set default smoothing parameters
  f_smooth.string <- paste(format(f_smooth), collapse='')

  if (f.smooth.string == paste(format(met_smooth_gam), collapse='')) {
    if (! 'k' %in% list(...)) {
      f_smooth <- function(...) f_smooth(..., k=5)
    }
  } else if (f.smooth.string == paste(format(met_smooth_loess), collapse='')) {
    if (! 'span' %in% list(...)) {
      f_smooth <- function(...) f_smooth(..., span=0.75)
    }
  }

  # Generating data frame, with a new "corrected" column
  d <- data.frame(time, concentration, corrected=concentration, compound)

  # Iteratively correcting systematic deviation
  for (i in 1:max.iter) {
    # Generating fit and calculating deviations
    d <- d %>%
           group_by(compound) %>%
           mutate(fit=f_smooth(time, corrected, ...),
                  deviation=(corrected-fit)/corrected)

    # Identifying the timepoint with largest deviation 
    deviations <- d %>%
                    group_by(time) %>%
                    summarize(med=abs(median(deviation, na.rm=TRUE))) %>%
                    arrange(desc(med))

    if (deviations$med[1] < min.deviation) break

    max_dev_time <- deviations$time[1]
   
    logic <- d$time == max_dev_time
    d[logic, ] <- within(d[logic, ], {
      corrected <- corrected - median(deviation)*corrected
    })
  }
  
  return(d$corrected)
}

