# Additional file 1 for BMC Systems Biology manuscript:
# Identification of internal standard addition error in 1H-NMR time-course data
# Stanislav Sokolenko, Marc G Aucoin
#
# This file contains functions used to perform the smooth correction on cell
# culture data.

library(dplyr)
library(mgcv)

#------------------------------------------------------------------------
# Defining smoothing function
f_smooth <- function(time, concentration) {
  # Smooths concentration data using a gam-5 model.
  #
  # Parameters
  # ----------
  # time: numeric
  #   x-values used for smoothing
  # concentration: numeric
  #   y-values used for smoothing
  
  d <- data.frame(time=time, concentration=concentration)
  model <- gam(concentration ~ s(time, bs='cr', k=5), data=d)
  fit <- as.vector(predict(model, d))

  return(fit)
}

#------------------------------------------------------------------------
correct_relative_deviation <- function(time, concentration, compound,
                                       max.iter=20, min.deviation=0.02) {
  # Performs iterative smoothing to identify and correct relative bias.
  #
  # Parameters
  # ----------
  # time: numeric
  #   x-values used for correction
  # concentration: numeric
  #   y-values used for correction
  # compound: character
  #   compound string corresponding to each x and y data point
  # max.iter: integer
  #   maximum number of iterative corrections to take (not needed with
  #   proper deviation threshold, but prevents infinite loops)
  # min.deviation: numeric
  #   minimum fractional deviation required to initiate correction

  # Generating data frame, with a new "corrected" column
  d <- data.frame(time, concentration, corrected=concentration, compound)

  # Iteratively correcting systematic deviation
  for (i in 1:max.iter) {
    # Generating fit and calculating deviations
    d <- d %>%
           group_by(compound) %>%
           mutate(fit=f_smooth(time, corrected),
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
  correction <- d$corrected-d$concentration
  
  return(correction)
}
