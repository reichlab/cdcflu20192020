#' Regional influenza incidence in the US (1997 - 2018)
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with 12,650 observations on weighted influenza-like illness 
#' measurements from all HHS regions, including the national level.
#' @source The cdcfluview R package. 
#' @docType data
#' @name flu_data
#' @usage data(flu_data)
NULL

#' Flu season "onset thresholds" from the US CDC.
#' 
#' @format A data.frame with 253 observations of the seasonal baseline threshold used for 
#' determining the "flu season onset" in each region of the US. Note that thresholds
#' for the 1997/1998 through 2006/2007 seasons have been imputed as the mean of later seasons.
#' @source \url{https://github.com/cdcepi/FluSight-forecasts/blob/master/wILI_Baseline.csv}
#' @docType data
#' @name flu_onset_baselines
#' @usage data(flu_onset_baselines)
NULL