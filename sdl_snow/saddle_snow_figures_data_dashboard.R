
# ============================================================================
# Niwot Ridge LTER Data Dashboard: Snow Depth Figures
# ============================================================================
# Initial code developed by Sarah C. Elmemdorf (2023)
# Code revised and updated for the data dashboard by Anne Marie Panetta (2026)

# Purpose:
#   Reads Saddle grid snow depth survey data and produces:
#
#     1. sdl_snow_anom.png   Annual May snow depth anomaly at Saddle
#
# Input:
#   data/saddsnow.dw.data.csv  -- Saddle grid snow depth surveys, one row
#                                  per grid point/date, with columns:
#                                  point_ID, date, mean_depth (some values
#                                  truncated at pole length, e.g. "210+")
#
# Output:
#   One PNG figure written to figures_dir (see file name above).
#
# TODO
#   - improve gap-filling for censored (e.g. >300cm) data
#   - check whether data sampling/coverage is temporally consistent for
#     May, or whether we want to pick a different month
#   - think about making better use of the full time series
#
# Original analysis: Sarah C. Elmendorf, 11 Jan 2023
# Revised and updated for the Data Dashboard: Anne Marie Panetta, July 2026
# ============================================================================

# -- SETUP --------------------------------------------------------------

rm(list = ls())

library(tidyverse)   # data wrangling + ggplot2
library(ggthemes)    # theme_hc()
library(lubridate)   # date handling
library(EDIutils)    # tools for interacting with EDI's data package API
library(ggpmisc)     # symmetric_limits() for scale_y_continuous

options(stringsAsFactors = FALSE)
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# Anchor working directory to this script's location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#source("../utility_functions/utility_functions_all.R")

data_dir <- "data"
figures_dir <- "figures"

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Set to TRUE to (re)download source data from EDI. Only needs to be run
# once, or when a newer package version is released.
download_data <- FALSE


# -- DOWNLOAD DATA FROM EDI ---------------------------------------------
# Note: if you have already downloaded SOME data, read_data_package_archive()
# will not overwrite existing files. Clear the /data directory first if you
# need to force a fresh download.
if (download_data) {
  # knb-lter-nwt.31 (Saddle grid snow depth), with current citation:
  #   Walker, D., J. Morse, and Niwot Ridge LTER. 2025. Snow depth data
  #   for Saddle grid, 1992 - ongoing. ver 22. Environmental Data
  #   Initiative. https://doi.org/10.6073/pasta/373278a850606e58d0cee389fbbad88b
  
  scope <- "knb-lter-nwt" # Niwot scope
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  for (id in c("31")) {
    # Ask EDI for the current version
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    packageID <- paste(scope, id, revision, sep = ".")
    
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID)) # confirm you're citing the version actually used
  }
  
  for (fname in list.files(data_dir, pattern = "knb-lter.*zip", full.names = TRUE)) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}

snowdata <- read_csv(file.path(data_dir, "saddsnow.dw.data.csv"))

# -- CLEAN DATA -----------------------------------------------------------

# It looks like a handful of measurements were skipped in a few years, but
# depth was 0 in that grid point both before and after, so we assume it
# was still 0 on the skipped date.
snowdata <- snowdata %>%
  mutate(
    mean_depth = ifelse(date == "2000-06-09" & point_ID == 301, 0, mean_depth),
    mean_depth = ifelse(date == "1998-01-09" & point_ID == 10, 0, mean_depth),
    mean_depth = ifelse(date == "1998-05-18" & point_ID == 401, 0, mean_depth)
  )

# Sampling seasons are often highly non-overlapping in date across years.
snowdata <- snowdata %>%
  mutate(
    doy = yday(date),
    season = ifelse(doy > 240, "fall", "spring"),
    year = lubridate::year(date)
  )


# seems like they skipped one meas in a few years, but before/after was 0
#in that grid so assume it's still 0
snowdata <- snowdata %>%
  mutate(
    mean_depth = ifelse(
      date == "2000-06-09" & point_ID == 301,
      0, mean_depth
    ),
    mean_depth = ifelse(
      date == "1998-01-09" & point_ID == 10,
      0, mean_depth
    ),
    mean_depth = ifelse(
      date == "1998-05-18" & point_ID == 401,
      0, mean_depth
    )
  )

# -- GAP-FILL CENSORED (210+/300+) DEPTH VALUES ----------------------------
# Some measurements are truncated at the pole length ("210+" or "300+"
# cm) rather than an exact depth. Fill these with the mean depth recorded
# at the same grid point within +/-15 days, restricted to other
# measurements that were also above the same threshold and greater than
# the next sample's depth (all such cases appear, on inspection, to be
# during the melt-down phase, so filling toward the higher/earlier end
# is appropriate here).
#
# Safe to ignore coercion warnings from as.numeric() on the "+"-suffixed
# strings below -- that's expected, not a data problem.

snowdata$fillval <- NA
for (i in 1:nrow(snowdata)) {
  if (snowdata$mean_depth[i] %in% c("210+", "300+") & !is.na(snowdata$mean_depth[i])) {
    mindoy <- snowdata$doy[i] - 15
    maxdoy <- snowdata$doy[i] + 15
    mypoint <- snowdata$point_ID[i]
    
    if (snowdata$mean_depth[i] == "210+" & !is.na(snowdata$mean_depth[i])) {
      fill_data <- snowdata %>%
        filter(doy > mindoy & doy < maxdoy & point_ID == snowdata$point_ID[i] & as.numeric(mean_depth) > 210)
    }
    if (snowdata$mean_depth[i] == "300+" & !is.na(snowdata$mean_depth[i])) {
      fill_data <- snowdata %>%
        filter(doy > mindoy & doy < maxdoy & point_ID == snowdata$point_ID[i] & as.numeric(mean_depth) > 300)
    }
    next_sample <- snowdata %>%
      filter(point_ID == snowdata$point_ID[i] & year == snowdata$year[i] & doy > snowdata$doy[i])
    next_sample <- next_sample %>%
      filter(doy == min(next_sample$doy))
    
    # All such cases, based on inspection, are in the melt-down phase when
    # this measurement was missed -- fill toward the higher/earlier bound.
    fill_data <- fill_data %>%
      filter(as.numeric(mean_depth) > as.numeric(next_sample$mean_depth))
    snowdata$fillval[i] <- mean(as.numeric(fill_data$mean_depth), na.rm = TRUE)
  }
}

# Spot-check the fill values against the original censored records.
probs <- snowdata[grepl("\\+", snowdata$mean_depth), ] %>%
  select(year, point_ID) %>%
  distinct()

probs <- probs %>%
  left_join(., snowdata) %>%
  select(year, doy, point_ID, mean_depth, fillval) %>%
  arrange(year, point_ID, doy)

# View(probs)  # looked ok on inspection


# -- GAP-FILL CENSORED (210+/300+) DEPTH VALUES ----------------------------
# Some measurements are truncated at the pole length ("210+" or "300+"
# cm) rather than an exact depth. Fill these with the mean depth recorded
# at the same grid point within +/-15 days, restricted to other
# measurements that were also above the same threshold and greater than
# the next sample's depth (all such cases appear, on inspection, to be
# during the melt-down phase, so filling toward the higher/earlier end
# is appropriate here).
#
# Safe to ignore coercion warnings from as.numeric() on the "+"-suffixed
# strings below -- that's expected, not a data problem.

snowdata$fillval <- NA
for (i in 1:nrow(snowdata)) {
  if (snowdata$mean_depth[i] %in% c("210+", "300+") & !is.na(snowdata$mean_depth[i])) {
    mindoy <- snowdata$doy[i] - 15
    maxdoy <- snowdata$doy[i] + 15
    mypoint <- snowdata$point_ID[i]
    
    if (snowdata$mean_depth[i] == "210+" & !is.na(snowdata$mean_depth[i])) {
      fill_data <- snowdata %>%
        filter(doy > mindoy & doy < maxdoy & point_ID == snowdata$point_ID[i] & as.numeric(mean_depth) > 210)
    }
    if (snowdata$mean_depth[i] == "300+" & !is.na(snowdata$mean_depth[i])) {
      fill_data <- snowdata %>%
        filter(doy > mindoy & doy < maxdoy & point_ID == snowdata$point_ID[i] & as.numeric(mean_depth) > 300)
    }
    next_sample <- snowdata %>%
      filter(point_ID == snowdata$point_ID[i] & year == snowdata$year[i] & doy > snowdata$doy[i])
    next_sample <- next_sample %>%
      filter(doy == min(next_sample$doy))
    
    # All such cases, based on inspection, are in the melt-down phase when
    # this measurement was missed -- fill toward the higher/earlier bound.
    fill_data <- fill_data %>%
      filter(as.numeric(mean_depth) > as.numeric(next_sample$mean_depth))
    snowdata$fillval[i] <- mean(as.numeric(fill_data$mean_depth), na.rm = TRUE)
  }
}

# Spot-check the fill values against the original censored records.
probs <- snowdata[grepl("\\+", snowdata$mean_depth), ] %>%
  select(year, point_ID) %>%
  distinct()

probs <- probs %>%
  left_join(., snowdata) %>%
  select(year, doy, point_ID, mean_depth, fillval) %>%
  arrange(year, point_ID, doy)

# View(probs)  # looked ok on inspection


# ============================================================================
# FIGURE 1: Annual May snow depth anomaly at Saddle
# ============================================================================
# Percent difference between each year's mean May snow depth and the
# long-term mean, colored by direction (above/below average).

snow_anom <- snowdata %>%
  dplyr::mutate(.,
                point_ID = as.integer(point_ID),
                date = lubridate::ymd(date),
                month = lubridate::month(date),
                year = lubridate::year(date),
                mean_depth = ifelse(grepl("\\+", mean_depth), fillval, mean_depth),
                mean_depth = as.numeric(as.character(mean_depth))
  ) %>%
  filter(month == 5) %>%
  group_by(year) %>%
  summarise(snowdepth = mean(mean_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    avg_depth = mean(snowdepth),
    depth_anom = snowdepth * 100 / avg_depth - 100,
    posneg = ifelse(depth_anom > 0, "pos", "neg") %>% factor(c("pos", "neg"))
  )

g1 <- ggplot(snow_anom, aes(x = year, y = depth_anom)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#4A9FE8", "#E85D8A")) +  # blue/pink, snowier vs less snowy
  labs(
    title = "May Snow Depth Anomalies",
    y = "Spring snow depth \n (% difference from long-term mean)",
    x = "Year"
  ) +
  scale_x_continuous(
    breaks = seq(
      5 * floor(min(snow_anom$year) / 5),
      5 * ceiling(max(snow_anom$year) / 5),
      by = 5
    )
  ) +
  scale_y_continuous(
    limits = c(-100, 100),
    sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(snowier %<->% less~snowy))
  ) +
  theme_hc() +
  theme(
    legend.position = "none"
  )
ggsave(g1,
       file = file.path(figures_dir, "sdl_snow_anom.png"),
       scale = 0.6, width = 10, height = 6
)

# Within-season anomaly (rather than a single May snapshot) was also
# explored, but sampling timing/coverage varies enough across years that
# gap-filling a within-season trajectory wasn't straightforward -- skipped
# for now.


# -- DERIVED METRICS ------------------------------------------------------
# The following blocks explore late-season gap-filling feasibility. They
# are not required for the figure above, but are kept here since later
# dashboard iterations may draw on them.

# Last spring survey that still had snow present.
no_snowfree <- snowdata %>%
  filter(season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(doy = max(doy)) %>%
  left_join(., snowdata) %>%
  filter(mean_depth != 0) %>%
  rename(last_snow_samp = doy) %>%
  arrange(year, point_ID, date)

# Last spring survey that was already snow-free -- these appear to fall
# late enough in the season (see doy note below) that no further
# gap-filling looks necessary for them.
lastdoy_is0 <- snowdata %>%
  filter(season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(doy = max(doy)) %>%
  left_join(., snowdata) %>%
  filter(mean_depth == 0) %>%
  rename(last_snow_samp = doy) %>%
  arrange(year, point_ID, date)

# Years/points where the last survey still had snow would need more
# complex gap-filling to resolve what happened afterward -- sampling
# appears to have simply stopped once the season got too snowy/late.
skip_stats <- no_snowfree %>%
  mutate(
    mean_depth = gsub("\\+", "", mean_depth),
    mean_depth = as.numeric(as.character(mean_depth))
  ) %>%
  group_by(year) %>%
  summarize(
    ct = dplyr::n(),
    mean_last_depth = mean(mean_depth),
    mean_last_doy = mean(last_snow_samp)
  )

# This late-season data might only be reliable up to about doy 153-157
# (already June 1), so it's not clear this adds much beyond the mean
# April/May snowpack already captured above.

