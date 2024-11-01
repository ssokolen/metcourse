# Functions for detecting systematic bias

#========================================================================>

#------------------------------------------------------------------------
#' B-spline fit
#'
#' Wraps bSpline function for use in bias detection and correction algorithm.
#'
#' @param time Time or sample number for metabolite time-course.
#' @param concentration Metabolite concentration.
#' @param ... Arguments to be passed into the bSpline() function,
#'            such as knots -- "The internal breakpoints that define the spline"
#'            and degree -- "Non-negative integer degree of the piecewise polynomial". If no additional arguments are provided, a reasonable
#'            default for is assumed.
#
#' @return A vector of fit concentrations.
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
#'
#' # Fitting
#' fit <- fit_b_spline(time, concentration,
#'                               knots = c(0.5), degree = 2)
#'
#' # Plotting
#' plot(time, concentration, type='p',
#'      xlab='Time post inoculation (hours)', ylab='Concentration (mM)')
#' lines(time, fit, type='l')
#' @export
fit_b_spline <- function(x, concentration, ...) {
  if (all(is.na(concentration))) return(concentration)

  d <- data.frame(x = x, concentration = concentration)

  f_basis <- function(d, knots = c(0.5), degree = 3) {
    B <- bSpline(d$x, knots = knots*max(d$x), degree = degree, intercept = TRUE)
  }

  if (length(list(...)) == 0) {
    B <- f_basis(d)
  } else {
    B <- f_basis(d, ...)
  }

  alpha <- solve(t(B) %*% B, t(B) %*% d$concentration)[,1]

  d$conc_pred <- B%*%alpha

  fit <- as.vector(d$conc_pred)

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
#' @param ... Arguments to be passed into \code{fit_b_spline}, (such as degree and number of knots).
#'
#' @return A dataframe corresponding to the time points with and without systematic error.
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
#' error <- detect_rel_bias(timecourse$time,
#'                                          timecourse$concentration,
#'                                          timecourse$metabolite)
#'
#'
#' @export
detect_rel_bias <- function(time, concentration, metabolite, min.deviation = NULL, ...) {

  d <- data.frame(time, concentration,metabolite)

  d <- d %>%
    group_by(metabolite) %>%
    arrange(time) %>%
    mutate(index = 1:length(time))

  extra_args <- list(...)

  # Passing in extra parameters into smoothing
  fit_b_spline_dply <- function(x, concentration) {
    args <- list(x = x, concentration = concentration)
    do.call(fit_b_spline, c(args, extra_args))
  }

  # Detecting systematic deviation
  # Generating fit and calculating deviations
  d <- d %>%
         group_by(metabolite) %>%
         mutate(fit = fit_b_spline_dply(time, concentration = concentration),
                deviation = (concentration - fit) / concentration)

  # Identifying the timepoint with large deviation
  deviations <- d %>%
                  group_by(time, index) %>%
                  summarize(med = abs(median(deviation, na.rm = TRUE))) %>%
                  arrange(desc(med)) %>%
                  mutate(rank = 1:n_distinct(time))

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

