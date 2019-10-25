## make a set of prospective prediction files for KDE/GAM model for 2017/2018 season
## Nicholas Reich
## October 2017: created
## October 2018: updated

library(cdcfluforecasts)
library(doMC)

FIRST_YEAR_OF_CURRENT_SEASON <- 2019
this_season <- paste0(FIRST_YEAR_OF_CURRENT_SEASON, "/", FIRST_YEAR_OF_CURRENT_SEASON+1)

data(flu_data)

n_sims <- 100000

## for 2018-2019 season, 
## first predictions due on 10/29/2018 (EW44 == SW14)  MMWRweek("2018-10-29")
## using data posted on 10/26/2016 that includes up to EW42 == SW12
## last predictions due on 5/13/2017 (EW20 == SW 42) MMWRweek("2019-05-13")
## last predictions use data through EW18 == SW40
## first_analysis_time_season_week could be set to 15, but padding at front end

season_weeks <- 10:43
region_strings <- c("National", paste("Region", 1:10))
fit_path <- "inst/estimation/region-kde/fits/"

registerDoMC(3)

## fit 2018/2019 models
foreach(reg=region_strings) %dopar% {
    
  ## fit models on training seasons, using only prospective data, not LOSO
  ## this function call saves a set of .rds files that contain the list defining a "KDE fit" 
  ## one fit for each (prospective season, region) pair
  
  # reg = region_strings[1]
  fit_region_kdes(flu_data, 
    region=reg,
    first_fit_year = FIRST_YEAR_OF_CURRENT_SEASON,
    last_fit_year = FIRST_YEAR_OF_CURRENT_SEASON,
    first_fit_week = 20, 
    path = fit_path)
}

## make entry files
foreach(season_week = season_weeks) %dopar% {
  ## season_week <- season_weeks[1]
  make_one_kde_prediction_file(save_path = "inst/submissions/region-kde/",
    fits_path = fit_path,
    season = this_season,
    season_week = season_week,
    n_sim = n_sims)
}

