# Functions for correcting systematic bias
#========================================================================>



#-----------------------------------------------------------------------
# Correction function

#------------------------------------------------------------------------
#' Correct systematic bias in metabolomic time-course data
#'
#' Identifies systematic deviations in metabolic data using a spline fit.
#' Timepoints that have a median relative deviation above a min.deviation value
#' are assumed to be influenced by a measurement or methodological bias as
#' compared to the overall trends metabolite concentrations.
#'
#' @param time A vector of times or sample numbers for metabolite time-courses.
#'             Note, there should be a time value for every \code{concentration} 
#'             for every \code{compound} e.g. time = c(1, 2, 3, 1, 2, 3), 
#'             concentration = c(20, 15, 10, 3, 6, 9), 
#'             compound = c('glc', 'glc', 'glc', 'lac', 'lac', 'lac'). 
#' @param concentration A vector of metabolite concentrations.
#' @param compound A vector of metabolite names that correspond to \code{time}
#'                   and \code{concentration}.
#' @param min.deviation Smallest median relative deviation to identify as a 
#'                      bias.
#' @param ... Arguments to be passed into \code{detect_rel_bias}, (such as knots
#'            and degree of the B-spline
#
#' @return A dataframe with the fit.
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
#' # Correcting
#' fit <- correct_rel_bias(timecourse$time, 
#'                                          timecourse$concentration,
#'                                          timecourse$metabolite)
#' timecourse$corrected <- fit$fit
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

correct_rel_bias <- function(time, compound, concentration, min.deviation = NULL, ...) { 
  #-----------------------------------------------------------------------------
  # Organize data
  
  d <- tibble(x = time, compound = compound, y = concentration)
  d <- d%>%
    arrange(compound, x)
  n.obs <- n_distinct(d$x)
  n.cmp <- n_distinct(d$compound)
  if (length(list(...)) == 0) {
    degree <- 3
    knots <- c(0.5)
  }
  # Assign a number to each compound 
  d$metabolite <- rep(1:n.cmp, each = n.obs) 
  
  #-----------------------------------------------------------------------------
  # Generate diagonal basis matrix 

  f_basis <- function(d) {
    B <- bSpline(d$x, knots = knots*max(d$x), degree = degree, intercept = TRUE)
  }

  f_basis_deriv <- function(d) {
    B <- dbs(d$x, knots = knots*max(d$x), degree = degree, intercept = TRUE)
  }
  
  bases <- d %>%
    split(d$metabolite) %>%
    map(f_basis)
  
  bases.deriv <- d %>%
    split(d$metabolite) %>%
    map(f_basis_deriv)
  
  B <- as.matrix(bdiagMat(bases))
  dB <- as.matrix(bdiagMat(bases.deriv))
  
  #-----------------------------------------------------------------------------
  # Calculate median deviation
  
  deviation <- detect_rel_bias(d$x, d$concentration, d$metabolite, min.deviation, ...)
  deviation <- deviation %>%
    filter(error == TRUE)

  #-----------------------------------------------------------------------------
  # Determine number of scaling terms that can be applied
  
  n.scal <- n.obs-(degree + 1 + length(knots)) 
  ind.scal <- deviation$index[1:min(n.scal, length(deviation$index))]

  # Checking if basis matrix is solvable
  count <- 0
  for (k in 1:min(n.scal, length(deviation$index))) {
    X <- B[1:n.obs, 1:(degree+1+length(knots))]
    X0 <- X
    actual <- svd(X0)
    actual <- actual$d
    min.actual <- min(actual)
    X <- X[-ind.scal, ]
    check <- svd(X)
    check <- check$d
    
    if (min(check) < 0.5*min.actual) { 
      ind.scal <- ind.scal[1:(length(ind.scal)-1)]
      count <- count + 1
    }
    else {
      break
    }
  }
  
  # Add terms (except the last one removed) back to check if solvable
  n.s <- length(ind.scal)
  c <- count - 1
  if (c > 0) {
    for (p in 1:c) {
      val <- deviation$index[(n.s+1+p)]
      ind.scal.store <- append(ind.scal, val)
      X <- B[1:n.obs, 1:(degree+1+length(knots))]
      X <- X[-ind.scal.store, ]
      check <- svd(X)
      check <- check$d
      if (min(check) > 0.5*min.actual) {
        ind.scal <- ind.scal.store
        count <- count - 1
      }
    }
  }
  # Ensure number of terms scaled doesn't exceed the degrees of freedom
  if (length(ind.scal) > n.scal) {
    ind.scal <- ind.scal[1:n.scal]
  } 
  #-----------------------------------------------------------------------------
  # Model fit
  
  # Parameters
  n <- nrow(d) 
  nb <- n.cmp*(degree + 1 + length(knots)) 
  d$time <- rep(1:n.obs, n.cmp)
  d$y[is.na(d$y)] <- 0
  
  input <- list(
    n = n, 
    nt = n.obs, 
    scal = scal, 
    sample = d$time, 
    nb = nb, 
    nc = n.cmp,
    compound = d$metabolite,
    B = B, 
    y = d$y
  )
  
  fit <- stan(file = "scaled_bias_function.stan", data = input , iter = 4000, control = list(max_treedepth = 10))

  #-----------------------------------------------------------------------------
  # Process output parameters
  
  par <- summary(fit)$summary[,"mean"]
  
  # B-spline parameters
  alpha.s <- par[1:nb]
  alpha.s <- unlist(alpha.s, use.names = FALSE)
  alpha.s <- matrix(alpha.s, nrow = nb)
  
  # Bias estimate
  bias.est <- par[(nb+n.obs+n.cmp+1):(nb+n.cmp+2*n.obs)]
  bias.est <- unlist(bias.est, use.names = FALSE)
  
  # Confidence intervals for bias estimate
  par.lower <- summary(fit)$summary[,"2.5%"]
  bias.est.2.5 <- par.lower[(nb+n.obs+n.cmp+1):(nb+n.cmp+2*n.obs)]
  bias.est.2.5 <- unlist(bias.est.2.5, use.names = FALSE)
  
  par.upper <- summary(fit)$summary[,"97.5%"]
  bias.est.97.5 <- par.upper[(nb+n.obs+n.cmp+1):(nb+n.cmp+2*n.obs)]
  bias.est.97.5 <- unlist(bias.est.97.5, use.names = FALSE)

  d$bias.est <- rep(bias.est, n.cmp)
  d$bias.est.2.5 <- rep(bias.est.2.5, n.cmp)
  d$bias.est.97.5 <- rep(bias.est.97.5, n.cmp)
  
  # Model fit
  d$fit <- B%*%alpha.s
  d$derivative <- dB%*%alpha.s

  # Output
  return(d)
}

