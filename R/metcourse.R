# Documentation for provided data-sets

#' Parameters used to simulate a set of realistic metabolic time-courses.
#'
#' Parameters designed to reflect commonly observed time-courses from a 
#' suspension culture of insect cells.
#'
#' @format A list with 24 parameters. See example of `simulate_timecourse`
#'         for more information.
"timecourse_param"

#' A set of realistic metabolic time-courses.
#'
#' A set of metabolic time-courses designed to mimic a suspension culture of 
#' insect cells. The parameters used to generate this data can be seen in 
#' `timecourse_param`.
#'
#' @format A data frame with 400 rows and 6 variables. The data simulates 1
#'         metabolomic experiment where 40 metabolites are tracked over 10
#'         timepoints.
#' \describe{
#'    \item{experiment}{Experiment number (a column of 1s in this case)}
#'    \item{metabolite}{Metabolites numbered 1 to 40}
#'    \item{sample}{Samples numbered 1 to 10}
#'    \item{concentration}{Metabolite concentrations (in mM)}
#'    \item{time}{Time post inoculation}
#'    \item{observed}{Metabolite concentrations with a known bias of 5% added 
#'                    at sample 4.}
#' }
"timecourse"


