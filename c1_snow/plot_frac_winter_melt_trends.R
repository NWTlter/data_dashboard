# Compare saddle snow to sumGDD
#SCE stopped on this about 7/11/2025
# setup
library(tidyverse)
library(snotelr)
library(trend) # for Mann-Kendall and Theil-Sen tests
library(gridExtra)
library(gtable)
library(grid)


# only need to download once
download_data_snotel <- FALSE
download_data_c1 <- FALSE

# functions -----------------------------------------------------------
# define function to perform
perform_trend_analysis <- function(data, value_col = "tot_ppt", year_col = "year") {
  # Remove any missing values
  clean_data <- data[!is.na(data[[value_col]]), ]
  
  if (nrow(clean_data) < 3) {
    return(data.frame(
      mk_tau = NA,
      mk_pvalue = NA,
      mk_trend = "insufficient data",
      ts_slope = NA,
      ts_intercept = NA
    ))
  }
  data <- data %>% arrange(year_col)
  # Mann-Kendall test
  mk_test <- mk.test(data[[value_col]])
  
  # Theil-Sen slope
  ts_test <- as.numeric(sens.slope(data[[value_col]])$estimates)
  
  # Calculate intercept
  intercepts <- data[[value_col]] - ts_test * data[[year_col]]
  sens_intercept <- median(intercepts)
  
  # Determine trend direction
  trend_direction <- case_when(
    mk_test$p.value > 0.05 ~ "no trend",
    mk_test$statistic > 0 ~ "increasing",
    mk_test$statistic < 0 ~ "decreasing",
    TRUE ~ "no trend"
  )
  
  return(data.frame(
    mk_tau = as.numeric(mk_test$statistic),
    mk_pvalue = as.numeric(mk_test$p.value),
    mk_trend = trend_direction,
    ts_slope = ts_test,
    ts_intercept = sens_intercept
  ))
}

# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data_c1) {
  # download the data from EDI
  # 185 is c1 temp
  
  scope <- "knb-lter-nwt" # Niwot scope
  
  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c("185"
  )) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    
    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")
    
    # download the data
    read_data_package_archive(packageID, path = "c1_d1_sdl_temp_ppt/data")
    print(read_data_package_citation(packageID))
  }
  
  # update the below so you remember to cite it correctly
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for C1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/0b4f8747ea72258e46681bb511c262f3. Accessed 2025-07-11."

  # unzip all zip files in the directory
  # overwrites the manifests but don't really need them.
  for (fname in list.files("c1_d1_sdl_temp_ppt/data",
                           pattern = "knb-lter.*zip",
                           full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = "c1_d1_sdl_temp_ppt/data/")
  }
}

if (download_data_snotel) {

# site_meta_data <- snotel_info()
# head(site_meta_data)

# but we know it's 838, so feed in manually
# snotel data?

# meta_data <- snotel_info() #if you need to sort out
# again other stations
# data starts Oct 1 1978
# but note there's no swe data til oct 1 1979

#often you have to run this twice I am not sure why
snotel_download(site_id = 838, path = 'c1_snow/data',
                                    internal = F)
}

# read and munge snow  --------------------------------------------------
snow_data_snotel<-read_csv('c1_snow/data/snotel_838.csv')


#check what dates avail
range(snow_data_snotel$date[!is.na(snow_data_snotel$snow_water_equivalent)])#so only water year 1980 - 2025ish good for snow at this point

# save data so as to not need to deal with finicky API again
saveRDS(snow_data_snotel, 'c1_d1_sdl_temp_ppt/data/snow_data_snotel.rds')

# calculate snow phenology,
#note this only runs if 365d in the year
# to the last/partial year will be cut off
snow_phenology_snotel <- snotel_phenology(snow_data_snotel)

#add the current year we are in if you wish
# max_2025<-max(snow_data_snotel$snow_water_equivalent[snow_data_snotel$date >= '2025-01-01'])
# max_2025_doy <- snow_data_snotel$date[snow_data_snotel$snow_water_equivalent == max_2025 &
#                         snow_data_snotel$date >= '2025-01-01'] %>%
#   lubridate::yday(.) %>% max(.)
# 
# snow_phenology_snotel <- snow_phenology_snotel %>%
#   bind_rows(., data.frame (year = 2025, max_swe = max_2025, max_swe_doy = max_2025_doy))




# Calculate frac melt pre max SWE -------------------------------------------
# Rationale from NCC paper
# We introduce a metric derived from daily SWE observations that complements
# the date of maximum SWE and provides additional hydrologic information.
# The fraction of cumulative annual snow water resources that has melted
# before a given date i, FMi, was computed for each station-year. This daily
# metric was computed in three steps for each of two dates: 1) the date of
# maximum SWE and 2) April 1st. First, daily melt (Figure S1; blue bars)
# was computed as the daily decrease in SWE, presented as positive values.
# Second, cumulative daily melt (Figure S1; red hashed line) was computed
# as the cumulative sum of daily melt from October 1st to August 1st.
# Third, the cumulative sum of daily melt was normalized by the total
# annual melt (on August 1st ) to estimate the (daily) fraction of cumulative
# annual melt, which was then sampled on the dates of maximum SWE (FMmax)
# and April 1st (FMApr1). The August 1st end date was chosen to avoid rare
# cases where early-season (i.e., September) snowfall could impact estimates
# of the total annual melt and to ensure that any late-lying snow (almost
# always gone by August at sensor locations) was recorded as annual melt.

frac_pre_swe <- snow_data_snotel %>%
  arrange(date) %>%
  # one missing value but it's july and 0 long before and after
  mutate(snow_water_equivalent = ifelse (date == '2022-07-16', 0, snow_water_equivalent),
         melt = snow_water_equivalent - lag(snow_water_equivalent, 1),
         melt = ifelse(melt > 0 & !is.na(melt), 0, melt),
         year = lubridate::year(date),
         month = lubridate::month(date),
         doy = lubridate::yday(date),
         water_year = ifelse(month >= 10, year + 1, year)) %>%
  inner_join(., snow_phenology_snotel) %>%
  filter(water_year >= 1980 & water_year < 2025) %>%
  mutate(bef = ifelse(doy < max_swe_doy | month >= 10, "bef", "after")) %>%
  #per keith's methods 0 out any melt in Aug, Sept
  # this doesn't actually matter here...
  mutate (melt = ifelse (month %in% c(8,9), 0, melt)) %>%
  group_by(bef, water_year) %>%
  summarize(melt_acc = -1 * sum(melt)) %>%
  pivot_wider(names_from = bef, values_from = melt_acc) %>%
  mutate(tot = after + bef,
         frac = bef / tot)


#calc fraction melted before april 1 instead?
prior_melt_apr_1 <- snow_data_snotel %>%
  arrange(date) %>%
  # one missing value but it's july and 0 long before and after
  mutate(snow_water_equivalent = ifelse (date == '2022-07-16', 0, snow_water_equivalent),
         melt = snow_water_equivalent - lag(snow_water_equivalent, 1),
         melt = ifelse(melt > 0 & !is.na(melt), 0, melt),
         year = lubridate::year(date),
         month = lubridate::month(date),
         doy = lubridate::yday(date),
         water_year = ifelse(month >= 10, year + 1, year)) %>%
  inner_join(., snow_phenology_snotel) %>%
  filter(water_year >= 1980 & water_year < 2025) %>%
  mutate(bef = ifelse(month <= 3 | month >= 10, "bef", "after")) %>%
  #per keith's methods 0 out any melt in Aug, Sept
  # this doesn't actually matter here...
  mutate (melt = ifelse (month %in% c(8,9), 0, melt)) %>%
  group_by(bef, water_year) %>%
  summarize(melt_acc = -1 * sum(melt)) %>%
  pivot_wider(names_from = bef, values_from = melt_acc) %>%
  mutate(tot = after + bef,
         frac = bef / tot)
  
#

# Calculate trends by site and season-------------------------------------------
trend_results_snowmelt_doy <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "last_snow_melt_doy", year_col = "year")
trend_results_max_swe_doy <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "max_swe_doy", year_col = "year")
trend_results_max_swe <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "max_swe", year_col = "year")
trend_results_max_frac_pre <- frac_pre_swe %>%
  perform_trend_analysis(., value_col = "frac", year_col = "water_year")
trend_results_prior_melt_apr_1 <- prior_melt_apr_1 %>%
  perform_trend_analysis(., value_col = "frac", year_col = "water_year")


# Create plot for a specific season with all sites
max_swe_doy <- create_seasonal_plot (seasonal_data = 
                           (snow_phenology_snotel %>%
                           mutate(season = 'annual',
                           site = 'c1')),
                                  trend_results = (trend_results_max_swe_doy %>%
                        mutate(site = 'c1', season = 'annual')),
                                  target_season = 'annual',
                                 y_col = "max_swe_doy", y_label = "max swe timing (doy)",
                                 plot_title = "Timing of max SWE", slope_units = "days")


snowmelt_doy <- create_seasonal_plot (seasonal_data = 
                                        (snow_phenology_snotel %>%
                                           mutate(season = 'annual',
                                                  site = 'c1')),
                                      trend_results = (trend_results_snowmelt_doy %>%
                                                         mutate(site = 'c1', season = 'annual')),
                                      target_season = 'annual',
                                      y_col = "last_snow_melt_doy", y_label = "last snow melt (doy)",
                                      plot_title = "Timing of last snowmelt", slope_units = "days")

max_swe <- create_seasonal_plot (seasonal_data = 
                                        (snow_phenology_snotel %>%
                                           mutate(season = 'annual',
                                                  site = 'c1')),
                                      trend_results = (trend_results_max_swe %>%
                                                         mutate(site = 'c1', season = 'annual')),
                                      target_season = 'annual',
                                      y_col = "max_swe", y_label = "max_swe (mm)",
                                      plot_title = "Max SWE", slope_units = "mm")

frac_pre <- create_seasonal_plot (seasonal_data = 
                                   (frac_pre_swe %>%
                                      mutate(season = 'annual',
                                             site = 'c1') %>%
                                      rename(year = water_year) %>%
                                      ungroup()),
                                 trend_results = (trend_results_max_frac_pre %>%
                                                    mutate(site = 'c1', season = 'annual')),
                                 target_season = 'annual',
                                 y_col = "frac", y_label = "fraction of melt \n before max swe",
                                 plot_title = "Fraction of winter melt",
                                 slope_units = "%")


# Extract legend from one plot
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

shared_legend <- get_legend(frac_pre)


jpeg('c1_snow/figures/combined_plot_snow.jpg', width = 12, height = 10, units = 'in', res = 300)
grid.arrange(max_swe_doy + theme(legend.position = "none"),
             snowmelt_doy + theme(legend.position = "none"),
             max_swe + theme(legend.position = "none"),
             frac_pre + theme(legend.position = "none"),
             shared_legend,
             ncol = 2, nrow = 3,
             layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
             heights = c(1, 1, 0.1)
)
dev.off()

