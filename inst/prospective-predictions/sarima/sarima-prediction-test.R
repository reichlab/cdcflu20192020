library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(MMWRweek)
library(cdcfluutils)
library(cdcflu20192020)
library(predx)
#library(FluSight)

### Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
analysis_time_year <- args[1]
analysis_time_week <- args[2]
seasonal_difference <- as.logical(args[3])
backfill_method <- args[4]

analysis_time_year <- "2018"
analysis_time_week <- "47"
seasonal_difference <- FALSE
backfill_method <- "none"
backfill_method <- "forecast_input"
#backfill_method <- "post-hoc"

method <- paste0("sarima_seasonal_difference_", seasonal_difference, "_backfill_", backfill_method)

submissions_save_path <- paste0("inst/prospective-predictions/sarima/prospective_predictions_", method)


analysis_time_epiweek <- paste0(
    analysis_time_year,
    sprintf("%02d", as.integer(analysis_time_week))) %>%
  as.integer()

## set up RNG
cdcflu20192020::get_rng_substream(
  method = method,
  year = analysis_time_year,
  week = analysis_time_week)

all_regions <- c("nat", paste0("hhs", 1:10))
data <- purrr::map_dfr(all_regions,
  function(region_str) {
    get_partially_revised_ilinet(
      region = region_str,
      epiweek = analysis_time_epiweek) %>%
    mutate(region = region_str)
  }) %>%
  mutate(
    region.type = ifelse(region == "nat", "National", "HHS Regions"),
    region = ifelse(region == "nat", "National", gsub("hhs", "Region ", region)),
    weighted_ili = wili,
    time = as.POSIXct(MMWRweek2Date(year, week))
  )

## Add time_index column: the number of days since some origin date (1970-1-1 in this case).
## The origin is arbitrary.
data$time_index <- as.integer(data$time -  as.POSIXct(ymd(paste("1970", "01", "01", sep = "-"))))

## Season column: for example, weeks of 2010 up through and including week 30 get season 2009/2010;
## weeks after week 30 get season 2010/2011
## I am not sure where I got that the season as defined as starting on MMWR week 30 from...
data$season <- ifelse(
  data$week <= 30,
  paste0(data$year - 1, "/", data$year),
  paste0(data$year, "/", data$year + 1)
)

## Season week column: week number within season
## weeks after week 30 get season_week = week - 30
## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
## This computation relies on the start_date function in package MMWRweek,
## which is not exported from that package's namespace!!!
data$season_week <- ifelse(
  data$week <= 30,
  data$week + MMWRweek(MMWRweek:::start_date(data$year) - 1)$MMWRweek - 30,
  data$week - 30
)

data <- as.data.frame(data)


analysis_time_season <- data$season[
  data$region == "National" &
  data$year == as.integer(analysis_time_year) &
  data$week == as.integer(analysis_time_week)
]


## Parameters used in simulating trajectories via sarima
simulate_trajectories_sarima_params <- list(
  fits_filepath = file.path(find.package("cdcflu20192020"),
    "estimation",
    "region-sarima",
    ifelse(seasonal_difference,
      "fits-seasonal-differencing",
      "fits-no-seasonal-differencing")),
  prediction_target_var = "weighted_ili",
  seasonal_difference = seasonal_difference,
  transformation = "box-cox",
  first_test_season = analysis_time_season,
  age = NA,
  regional_switch = "NatRegState"
)

weeks_in_first_season_year <-
  get_num_MMWR_weeks_in_first_season_year(analysis_time_season)



all_regions <- c(paste0("Region ", 1:10), "National")
res_predx <- purrr::map_dfr(all_regions, function(region) {
  get_submission_one_region_via_trajectory_simulation(
    data = data,
    analysis_time_season = analysis_time_season,
    first_analysis_time_season_week = 10, # == week 40 of year
    last_analysis_time_season_week = weeks_in_first_season_year - 11, # analysis for 33-week season
    region = region,
    prediction_target_var = "weighted_ili",
    incidence_bins = data.frame(
      lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
      upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
    n_trajectory_sims = 100,
    simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
    simulate_trajectories_params = simulate_trajectories_sarima_params,
    backfill_adjust = backfill_method) #, "forecast input", "post-hoc"
}) %>%
  mutate(
    location = ifelse(
      location == "National",
      "US National",
      paste0("HHS ", location)
    )
  )


predx_res_file <- file.path(submissions_save_path,
                          "predx",
                          paste0(
                            "EW", analysis_time_week,
                            "-", analysis_time_year,
                            "-ReichLab_", method,
                            ".rds"))
saveRDS(res_predx, predx_res_file)

res_csv <- res_predx %>%
  dplyr::filter(predx_class %in% c("BinCat", "BinLwr", "Point")) %>%
  dplyr::mutate(
    team = "Kernel of Truth",
    mmwr_week = as.character(analysis_time_week),
    submission_date = Sys.Date(),
    unit = ifelse(target %in% c("Season onset", "Season peak week"), "week", "percent")
  ) %>%
  predx::export_flusight_csv()

csv_res_file <- file.path(submissions_save_path,
  "csv",
  paste0(
    "EW", analysis_time_week,
    "-", analysis_time_year,
    "-ReichLab_", method,
    ".csv"))
write.csv(res_csv,
  file = csv_res_file,
  row.names = FALSE)

#file_valid <- FluSight::verify_entry_file(csv_res_file)
#
#if(!file_valid) {
#  stop("invalid predictions file!")
#}
