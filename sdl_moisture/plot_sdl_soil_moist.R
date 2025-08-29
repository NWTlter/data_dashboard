# Script to plot sdl (nwt package id 405) soil moisture
# temp patterns

library(tidyverse)
library(ggthemes)
library(lemon)

rm(list = ls())

# only need to download once
download_data <- FALSE

# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data) {
  # download the data from EDI
  # 405 is sdl temp
  scope <- "knb-lter-nwt" # Niwot scope
  
  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c(
    "405"
  )) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    
    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")
    
    # download the data
    read_data_package_archive(packageID, path = "sdl_moisture/data")
    print(read_data_package_citation(packageID))
  }
  
  # update the below so you remember to cite it correctly
  # "Morse, J. and M. Losleben. 2025. Climate data for saddle data loggers (CR23X and CR1000), 2000 - ongoing, daily. ver 10. Environmental Data Initiative. https://doi.org/10.6073/pasta/b01aea637f7608f0a1b2895ae474d571. Accessed 2025-08-11."
  # overwrites the manifests but don't really need them.
  for (fname in list.files("sdl_moisture/data",
                           pattern = "knb-lter.*zip",
                           full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = "sdl_moisture/data/")
  }
}

# read and munge ppt and temp --------------------------------------------------

df <- read.csv("sdl_moisture/data/sdlcr23x-cr1000.daily.ml.data.csv",
  na.strings = "NaN"
)

# remove potentially bad data from the soil temp & moist
# note flagging not nec done in all yrs
# so may want to redo some qc on the earlier yrs
# at some point

df <- df %>%
  mutate(
    date = lubridate::ymd(date),
    soiltemp_5cm_max = ifelse(grepl("q", flag_soiltemp_5cm_max), NA, soiltemp_5cm_max),
    soiltemp_5cm_min = ifelse(grepl("q", flag_soiltemp_5cm_min), NA, soiltemp_5cm_min),
    soilmoist_5cm_avg = ifelse(grepl("q", flag_soilmoist_5cm_avg), NA, soilmoist_5cm_avg),
    month = lubridate::month(date)
  )


# linearly interpolate up to 10d gaps for moisture
df <- df %>%
  arrange(date) %>%
  mutate(
    soiltemp_5cm_min = zoo::na.approx(soiltemp_5cm_min, maxgap = 10, rule = 2),
    soilmoist_5cm_avg = zoo::na.approx(soilmoist_5cm_avg, maxgap = 10, rule = 2)
  ) %>%
  # for gs chai calcs if we want to do those
  mutate(tmin = (soiltemp_5cm_min - -2) / (5 - -2)) %>%
  mutate(tmin = ifelse(tmin < 0, 0, tmin)) %>%
  mutate(tmin = ifelse(tmin > 1, 1, tmin)) %>%
  mutate(mmin = (soilmoist_5cm_avg - 0.05) / (.15 - 0.05)) %>%
  mutate(mmin = ifelse(mmin < 0, 0, mmin)) %>%
  mutate(mmin = ifelse(mmin > 1, 1, mmin))

summer_moist <- df %>%
  filter(month %in% c(6, 7, 8)) %>%
  group_by(year) %>%
  summarise(
    summer_soil_moist = mean(soilmoist_5cm_avg, na.rm = TRUE),
    ct = dplyr::n(), .groups = NULL
  ) %>%
  filter(ct > 80) %>% # make sure it's pretty complete after interpo
  mutate(
    avgmoist = mean(summer_soil_moist, na.rm = TRUE),
    anom_moist = (summer_soil_moist * 100 / avgmoist) - 100,
    posneg = ifelse(anom_moist > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )
# plot
g1 <- ggplot(summer_moist, aes(x = year, y = anom_moist)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#034e7b", "#99000d")) +
  # labs(y = expression(atop(paste("Distribution of node ages, ",
  #                                 italic(gamma)))))
  labs(y = "Soil moisture anomaly (%) \n Jun-Aug", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(wetter %<->% drier))) +
  theme_hc() +
  theme(legend.position = "none")

ggsave(g1,
  file = "sdl_moisture/figures/soil_moist_anom.png",
  scale = 0.5, width = 8, height = 6
)


summer_moist_last_yr <- df %>%
  mutate(yday = lubridate::yday(date)) %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  filter(year == max(df$year))

byday <- df %>%
  mutate(yday = lubridate::yday(date)) %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  filter(year != max(df$year)) %>%
  group_by(yday) %>%
  summarise(
    q10 = quantile(soilmoist_5cm_avg, 0.1, na.rm = TRUE),
    q50 = quantile(soilmoist_5cm_avg, 0.5, na.rm = TRUE),
    q90 = quantile(soilmoist_5cm_avg, 0.9, na.rm = TRUE),
    .groups = 'drop'
  )

fill_legend <- paste0("10th-90th percentile \n","(",
       min(df$year),"-",
       max(df$year)-1, ")")

g2 <- ggplot(summer_moist_last_yr %>%
  ungroup(), aes(x = yday)) +
  geom_line(aes(y=soilmoist_5cm_avg,
                color = as.character(year)))+
  geom_ribbon(data = byday, 
              aes(x = yday, ymin = q10, ymax = q90, 
                  fill = fill_legend), 
              alpha = 0.3)+
  scale_x_continuous(
    breaks = c(121, 152, 182, 213, 244),  # May 1, June 1, July 1, Aug 1, Sep 1 (approx)
    labels = c("May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1")
  ) +
  scale_color_manual(
    name = NULL,
    values = setNames("black", as.character(unique(summer_moist_last_yr$year)))
  )+
  scale_fill_manual(
    name = NULL,
    values = setNames("pink", fill_legend)
  )+
  ylab('Soil moisture (%)')+
  xlab(NULL)+
  theme_hc()

  
ggsave(g2,
  file = "sdl_moisture/figures/sdl_soil_moisture_anom.png",
  scale = 0.5, width = 8, height = 6
)

