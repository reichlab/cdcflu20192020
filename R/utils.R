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
