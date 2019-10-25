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
