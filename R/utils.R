#' Download and preprocess the latest CDC flu data, both national and regional
#'
#' @param latest_year year through which data should be downloaded, defaults to current year
#'
#' @return data frame with latest flu data, preprocessed
#' @export
download_and_preprocess_flu_data <- function(latest_year = as.numeric(format(Sys.Date(), "%Y"))) {  
  require(cdcfluview)
  require(lubridate)
  require(dplyr)
  require(MMWRweek)
  
  regionflu <- ilinet(region="hhs", years= 1997:latest_year)
  regionflu$region <- as.character(regionflu$region)
  
  usflu <- ilinet(region="national", years= 1997:latest_year)
  
  flu_data <- bind_rows(regionflu, usflu)
  
  flu_data <- transmute(flu_data,
    region_type = region_type,
    region = as.factor(region),
    year = year,
    week = week,
    time = as.POSIXct(MMWRweek2Date(year, week)),
    weighted_ili = weighted_ili)
  
  ## set zeroes to NAs
  flu_data[which(flu_data$weighted_ili==0),"weighted_ili"] <- NA
  
  ## Add time_index column: the number of days since some origin date
  ## (1970-1-1 in this case).  The origin is arbitrary.
  flu_data$time_index <- as.integer(lubridate::date(flu_data$time) -  ymd("1970-01-01"))
  
  ## Season column: for example, weeks of 2010 up through and including week 30
  ## get season 2009/2010; weeks after week 30 get season 2010/2011
  ## Official CDC flu season for the purposes of prediction runs from week 40 of
  ## one year to week 20 of the next; the season start week we define here is the
  ## mid-point of the "off-season"
  flu_data$season <- ifelse(
    flu_data$week <= 30,
    paste0(flu_data$year - 1, "/", flu_data$year),
    paste0(flu_data$year, "/", flu_data$year + 1)
  )
  
  ## Season week column: week number within season
  ## weeks after week 30 get season_week = week - 30
  ## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
  ## This computation relies on the start_date function in package MMWRweek,
  ## which is not exported from that package's namespace!!!
  flu_data$season_week <- ifelse(
    flu_data$week <= 30,
    flu_data$week + MMWRweek(MMWRweek:::start_date(flu_data$year) - 1)$MMWRweek - 30,
    flu_data$week - 30
  )
  
  flu_data <- as.data.frame(flu_data)
  
  return(flu_data)
}

#' Download and preprocess the latest CDC flu data, state-level
#'
#' @param latest_year year through which data should be downloaded, defaults to current year
#'
#' @return data frame with latest state-level flu data, preprocessed
#' @export
download_and_preprocess_state_flu_data <- function(latest_year = as.numeric(format(Sys.Date(), "%Y"))) {
  
  require(cdcfluview)
  require(MMWRweek)
  require(dplyr)
  require(lubridate)
  
  flu_data_raw <- ilinet(region="state", years=1997:latest_year)
  
  flu_data <- mutate(flu_data_raw, time = as.POSIXct(MMWRweek2Date(year, week)))
  
  ## set rows with denominator zeroes to NAs
  flu_data[which(flu_data$total_patients==0),"weighted_ili"] <- NA
  
  ## Add time_index column: the number of days since some origin date
  ## (1970-1-1 in this case).  The origin is arbitrary.
  flu_data$time_index <- as.integer(lubridate::date(flu_data$time) -  ymd("1970-01-01"))
  
  ## Season column: for example, weeks of 2010 up through and including week 30
  ## get season 2009/2010; weeks after week 30 get season 2010/2011
  ## Official CDC flu season for the purposes of prediction runs from week 40 of
  ## one year to week 20 of the next; the season start week we define here is the
  ## mid-point of the "off-season"
  flu_data$season <- ifelse(
    flu_data$week <= 30,
    paste0(flu_data$year - 1, "/", flu_data$year),
    paste0(flu_data$year, "/", flu_data$year + 1)
  )
  
  ## Season week column: week number within season
  ## weeks after week 30 get season_week = week - 30
  ## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
  ## This computation relies on the start_date function in package MMWRweek,
  ## which is not exported from that package's namespace!!!
  flu_data$season_week <- ifelse(
    flu_data$week <= 30,
    flu_data$week + MMWRweek(MMWRweek:::start_date(flu_data$year) - 1)$MMWRweek - 30,
    flu_data$week - 30
  )
  
  state_flu <- as.data.frame(flu_data)
  
  return(state_flu)
}


#' Get the initial rng substream for an rstream object.  Should be called by all
#' prediction methods before doing any random number generation.
#'
#' This function DOES NOT set RNG to use rstream with the returned object.
#' Because of strange behavior in the rstream package, this function can be
#' called only once in a given R session.
#'
#' @param seed integer seed for rng; the default was randomly generated
#'
#' @return object of class "rstream.mrg32k3a".  The object has been packed via
#'   rstream.packed.
#'
#' @export
get_initial_rng_substream <- function(
  seed = 9029979) {
  require("rstream")
  
  set.seed(seed)
  rngstream <- new("rstream.mrg32k3a", seed = sample(1:100000, 6, rep = FALSE))
  
  ## pack rngstream object and return (invisibly) in case methods want to use
  rstream.packed(rngstream) <- TRUE
  return(rngstream)
}


#' Get the rng substream for an rstream object corresponding to the combination
#' of prediction method, region, and season left out.  Should be called by all
#' prediction methods before doing any random number generation.
#'
#' Importantly, by default this function has the side effect of setting RNG to
#' use rstream with the returned object.  This behavior is determined by the
#' set_rng argument.  This means the caller doesn't have to worry about doing
#' anything unless (a) it wants to use more than 1 substream or (b) it is going
#' to parallelize or do any RNG in a different R session.  Because of strange
#' behavior in the rstream package, this function can be called at most once
#' without the rngstream argument in a given R session.
#'
#' @param rngstream (optional) object of class "rstream.mrg32k3a" which will be
#'   advanced from its current state.
#' @param seed integer seed for rng; the default was randomly generated
#' @param method character string specifying prediction method
#'   currently one of "sarima", "kcde", or "kde"
#' @param year character string specifying year, format "1998"
#' @param week character string specifying epidemic week, in format "02"
#' @param set_rng boolean should rng be set to use rstream with the returned
#'   rngstream object?
#'
#' @return invisible object of class "rstream.mrg32k3a", but advanced to the
#' first substream reserved for the given combination of prediction method,
#' region, and season.  The object has been packed via rstream.packed.
#'
#' @export
get_rng_substream <- function(
  rngstream,
  seed = 9029979,
  method,
  year,
  week,
  set_rng = TRUE) {
  require("rstream")
  
  ## Create a data frame with combinations of method, year and week,
  ## number of substreams used for each such combination.
  ## We can add more methods later without causing any problems by appending
  ## new method names to the END of the "method" vector below.
  ## Adding new years or weeks must be done by adding a new set of rows
  ## to the bottom of the substreams_used data frame (e.g. via bind_rows).
  year_week_combos <- expand.grid(
    year = as.character(2010:2019),
    week = sprintf("%02d", c(1:20, 40:52)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(epiweek = as.integer(paste0(year, week))) %>%
    filter(epiweek >= 201040 &
             epiweek <= 201920) %>%
    rbind(
      data.frame(year = "2014",
                 week = "53",
                 epiweek = 201453,
                 stringsAsFactors = FALSE)
    ) %>%
    arrange(epiweek)
  
  substreams_used <- expand.grid(
    epiweek = year_week_combos$epiweek,
    method = c(
      paste0("sarima_seasonal_difference_TRUE_",
        c("backfill_none", "backfill_forecast_input", "backfill_post-hoc")),
      paste0("sarima_seasonal_difference_FALSE_",
        c("backfill_none", "backfill_forecast_input", "backfill_post-hoc")),
      paste0("kcde",
        c("backfill_none", "backfill_forecast_input", "backfill_post-hoc")),
      "kde"),
    stringsAsFactors = FALSE
  )
  substreams_used$num_substreams <- 1
  ## if any future method uses more than 1 substream, set that here
  
  ## substream index for the specified method, region, and season
  ind <- which(
    substreams_used$epiweek == paste0(year, week) &
      substreams_used$method == method)
  
  if(!identical(length(ind), 1L)) {
    stop("Invalid year, week, and/or method.")
  }
  
  ## Create Rstream object and advance past all substreams used by previous
  ## methods/regions/seasons
  if(missing(rngstream)) {
    set.seed(seed)
    rngstream <- new("rstream.mrg32k3a", seed = sample(1:100000, 6, rep = FALSE))
  } else {
    rstream.packed(rngstream) <- FALSE
  }
  
  advance_count <- sum(substreams_used$num_substreams[seq_len(ind - 1)])
  for(i in seq_len(advance_count)) {
    rstream.nextsubstream(rngstream)
  }
  
  ## set to use rstream package for RNG with rngstream object
  if(set_rng) {
    rstream.RNG(rngstream)
  }
  
  ## pack rngstream object and return (invisibly) in case methods want to use
  rstream.packed(rngstream) <- TRUE
  invisible(rngstream)
}

#' Compute season onset, peak week, and peak incidence
#'
#' @param data a data frame containing at minimum columns named season,
#'   season_week and a column with some sort of incidence measure
#' @param season the season to look at
#' @param first_CDC_season_week the first week of the season to use for
#'   calculating onset and peak
#' @param last_CDC_season_week the last week of the season to use for
#'   calculating onset and peak
#' @param onset_baseline numeric baseline value for determining season onset
#' @param incidence_var a character string naming the variable in the data
#'   argument containing a measure of incidence, or an integer index
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#'
#' @return a list with four entries:
#'   1) observed_onset_week, either an integer between first_CDC_season_week
#'     and last_CDC_season_week (inclusive), or "none"
#'   2) observed_peak_week, an integer between first_CDC_season_week and
#'     last_CDC_season_week (inclusive)
#'   3) observed_peak_inc, a numeric with the maximum value of the specified
#'     incidence measure between first_CDC_season_week and last_CDC_season_week
#'   4) observed_peak_inc_bin, character name of incidence bin for peak incidence
#'
#' @export
get_observed_seasonal_quantities <- function(
  data,
  season,
  first_CDC_season_week = 10,
  last_CDC_season_week = 42,
  onset_baseline,
  incidence_var,
  incidence_bins,
  incidence_bin_names
) {
  first_season_ind <- min(which(data$season == season))
  last_season_ind <- max(which(data$season == season))
  
  obs_inc_in_season_leading_trailing_nas <-
    data[seq(from = first_season_ind, to = last_season_ind),
      incidence_var]
  
  ## pad so that we start at season week 1
  if(data$season_week[first_season_ind] != 1) {
    obs_inc_in_season_leading_trailing_nas <- c(
      rep(NA, data$season_week[first_season_ind] - 1),
      obs_inc_in_season_leading_trailing_nas)
  }
  
  ## set values before first analysis time season week or after last
  ## analysis time season week to NA
  ## these are outside of the bounds of the season the CDC wants to look at
  obs_inc_in_season_leading_trailing_nas[
    seq_len(first_CDC_season_week - 1)] <- NA
  if(length(obs_inc_in_season_leading_trailing_nas) >
      last_CDC_season_week) {
    obs_inc_in_season_leading_trailing_nas[
      seq(from = last_CDC_season_week + 1,
        to = length(obs_inc_in_season_leading_trailing_nas))] <- NA
  }
  
  observed_peak_inc <- max(
    obs_inc_in_season_leading_trailing_nas,
    na.rm = TRUE)
  observed_peak_inc_bin <- get_inc_bin(observed_peak_inc, return_character = TRUE)
  
  ## peak week timing is based on rounded values
  round_to_.1 <- function(inc_val) {
    if(is.na(inc_val)) {
      return(inc_val)
    } else {
      floor_val <- floor(inc_val * 10) / 10
      if(inc_val >= floor_val + 0.05) {
        return(floor_val + 0.1)
      } else {
        return(floor_val)
      }
    }
  }
  
  rounded_observed_peak_inc <- round_to_.1(observed_peak_inc)
  rounded_obs_inc_in_season <- sapply(obs_inc_in_season_leading_trailing_nas,
    round_to_.1
  )
  
  observed_peak_week <-
    which(rounded_obs_inc_in_season == as.numeric(rounded_observed_peak_inc))
  
  weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(season)
  observed_onset_week <- get_onset_week(
    incidence_trajectory = rounded_obs_inc_in_season,
    #    incidence_trajectory = obs_inc_in_season_leading_trailing_nas, # used in stable method
    baseline = onset_baseline,
    onset_length = 3L,
    first_season_week = 31,
    weeks_in_first_season_year = weeks_in_first_season_year
  )
  
  return(list(observed_onset_week = observed_onset_week,
    observed_peak_week = observed_peak_week,
    observed_peak_inc = observed_peak_inc,
    observed_peak_inc_bin = observed_peak_inc_bin
  ))
}


#' return the bin name for a given incidence
#'
#' @param inc numeric incidence level
#' @param return_character logical: if true, return type is character (bin name)
#'   if false, return type is numeric representation of bin
#'
#' @return vector giving the bin name of the input incidence.
#'
#' @details assumes max inc bin is 13 and bins are 0.1 in size.
#'
#' @export
get_inc_bin <- function(inc,max=13,
  return_character = TRUE) {
  inc <- round(inc, 1)
  bin_numeric <- ifelse(inc < max,
    floor(inc*10)/10, ## floors to 1st decimal place
    max)
  if(return_character) {
    return(as.character(bin_numeric))
  } else {
    return(bin_numeric)
  }
}

#' return integer that's either 52 or 53: number of weeks in the first year of
#' a given season.
#'
#' @param season season in the format "2014/2015"
#'
#' @details requires MMWRweek package
#'
#' @export
get_num_MMWR_weeks_in_first_season_year <- function(season) {
  return(get_num_MMWR_weeks_in_year(substr(season, 1, 4)))
}

#' return integer that's either 52 or 53: number of MMWR weeks in the given year
#'
#' @param year year in the format "2014" -- can be character or numeric
#'
#' @details requires non-exported function start_date from MMWRweek package
#'
#' @export
get_num_MMWR_weeks_in_year <- function(year) {
  require(MMWRweek)
  year <- as.numeric(year)
  return(MMWRweek::MMWRweek(MMWRweek:::start_date(year + 1) - 1)$MMWRweek)
}

#' Utility function to compute onset week based on a trajectory of incidence values
#'
#' @param incidence_trajectory a numeric vector with incidence for each time
#'   point in a season
#' @param baseline the threshold that incidence must cross to count as an onset
#' @param onset_length number of consecutive time points that incidence must
#'   exceed the baseline threshold in order to count as the season onset
#' @param first_season_week number of weeks in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' @param weeks_in_first_season_year How many MMWR weeks are in the first year
#'   of the season?  For example, in the 2000/2001 season, the first year is
#'   2000.  There were 52 MMWR weeks in 2000.
#'
#' @return the smallest index i such that every entry of
#'   incidence_trajectory[seq(from = i, length = onset_length)]
#'   is >= baseline, if such an index exists.
#'   Otherwise, the character vector "none"
#'
#' @export
get_onset_week <- function(incidence_trajectory,
  baseline,
  onset_length,
  first_season_week = 31,
  weeks_in_first_season_year) {
  
  exceeded_threshold <- sapply(
    seq_len(length(incidence_trajectory) - onset_length),
    function(start_ind) {
      above_baseline <- incidence_trajectory[seq(from = start_ind, length = onset_length)] >= baseline
      length(above_baseline)>0 &&
        all(above_baseline) &&
        !all(is.na(incidence_trajectory))
    }
  )
  
  if(any(exceeded_threshold, na.rm = TRUE)) {
    season_week <- min(which(exceeded_threshold))
    
    return(season_week)
  } else {
    return("none")
  }
}


#' Get the onset baseline for a combination of region and season
#'
#' @param region a string, either "National", "Region k", or "Regionk" where
#'   k in {1, ..., 10}
#' @param season a string, in the format "2015/2016"
#'
#' @return baseline value for determining season onset
#'
#' @export
get_onset_baseline <- function(region, season = "2015/2016") {
  ## pick baseline
  ## assumes region is either "National" or "Region k" format
  reg_string <- ifelse(region=="National", "National", gsub(" ", "", region))
  idx <- which(flu_onset_baselines$region==reg_string &
      flu_onset_baselines$season==season)
  reg_baseline <- flu_onset_baselines[idx, "baseline"]
  
  return(reg_baseline)
}

#' Calcluation of median value from binned probability distribution
#'
#' @param probs vector of named probabilities
#'
#' @return a numeric value
#'
#' @export
calc_median_from_binned_probs <- function(probs) {
  ## could do something more intelligent for "none" bin in onset - currently assuming it is all ordered
  cumprob <- cumsum(probs)
  median_idx <- min(which(cumprob>=0.5))
  as.numeric(names(probs)[median_idx])
}


