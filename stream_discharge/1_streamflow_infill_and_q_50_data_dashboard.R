# ============================================================================
# Niwot Ridge LTER Data Dashboard: Streamflow Gap-Filling and Export
# ============================================================================
#
# Purpose:
#   Downloads (optionally) and gap-fills daily discharge data for four
#   Niwot Ridge stream gauges -- Albion (alb), Green Lake 4 (gl4),
#   Martinelli (mar), and Saddle (sdl) -- and exports a combined,
#   gap-filled dataset for use in downstream plotting scripts.
#
#   Gaps are filled using a combination of:
#     - Known dry-channel periods (from site notes), set to zero flow
#     - Multivariate EM imputation (mnimput), using each site's other
#       years and/or correlated nearby sites as predictors
#     - Short (<=3 day) linear interpolation for any remaining small gaps
#
# Input:
#   EDI data packages (downloaded if download_data = TRUE; otherwise
#   expects the raw CSVs to already exist in data_dir -- see the
#   read-in section below for expected filenames)
#
# Output:
#   data/spctl_<sp_control>.csv -- combined, gap-filled daily discharge
#   data for all four sites, with columns: local_site, date, discharge,
#   yday, year, wyear, is_infilled
#
#   This file is the required input for streamflow_data_dashboard.R.
#
# Original analysis: Sarah C. Elmendorf, 20 Jan 2021
# Data Dashboard update: Anne Marie Panetta, 2026
#   (2026 update: reads directly from the latest EDI package version;
#   see citations below)
# ============================================================================

# todo (as of 2026 update):
#   - Consider retaining records currently excluded because they are
#     flagged "estimate based on partial records only"
#   - Consider weighting downstream analyses by amount of infilled data,
#     or excluding years that are mostly infilled


# -- SETUP ------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(mtsdi) # for imputing

# Anchor working directory to this script's location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Spline degrees-of-freedom control used by mnimput. Values of 3, 5, 7,
# and 9 were compared visually during development; 7 was chosen, but the
# final infilled values are not very sensitive to this choice.
sp_control <- 7

# If TRUE, generates in-session diagnostic plots at each infilling step
# so you can visually check how well the imputation performed. These are
# not saved to disk; set to TRUE when reviewing/debugging the pipeline.
plot_in_code <- FALSE

# Set to TRUE to (re)download source data from EDI. Only needs to be run
# once, or when a newer package version is released.
download_data <- FALSE

data_dir <- "data"
figures_dir <- "figures"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# -- DOWNLOAD DATA FROM EDI ---------------------------------------------
# Note: if you have already downloaded SOME data, read_data_package_archive()
# will not overwrite existing files. Clear the /data directory first if you
# need to force a fresh download.

if (download_data) {
  # Discharge data packages used, with current citations:

  # knb-lter-nwt.102 (albion)
    #Caine, N., J. Morse, and Niwot Ridge LTER. 2025. 
    #Streamflow data for Albion camp, 1981 - ongoing. ver 20. 
    #Environmental Data Initiative. 
    #https://doi.org/10.6073/pasta/00341116ab5c8eac60e641cb1b5c3468
  # knb-lter-nwt.105 (gl4)
    #Caine, N., J. Morse, and Niwot Ridge LTER. 2025. 
    #Streamflow for Green Lake 4, 1981 - ongoing. ver 20. 
    #Environmental Data Initiative. 
    #https://doi.org/10.6073/pasta/c3fa75cfe47c10fb61d866fe4d75a93a.
  # knb-lter-nwt.111 (martinelli)
    #Caine, T., J. Morse, and Niwot Ridge LTER. 2026. 
    #Streamflow for Martinelli basin, 1982 - ongoing. ver 17. 
    #Environmental Data Initiative. 
    #https://doi.org/10.6073/pasta/44b350d46f371082b2cc98490cf36959
  # knb-lter-nwt.74 (saddle)
    #Caine, T., J. Morse, and Niwot Ridge LTER. 2026. 
    #Streamflow data for Saddle stream, 1999 - ongoing. ver 10. 
    #Environmental Data Initiative. 
    #https://doi.org/10.6073/pasta/c699680482443efff3ad9f30c2c962a1
  
  # Related packages not used in this pipeline, but available on EDI if
  # a future revision wants additional sites: 163 (Green Lake 5 rock
  # glacier), 162 (watershed flume), 213 (Soddie), 109 (Green Lake 5)

  
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  scope <- "knb-lter-nwt" # Niwot scope
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  for (id in c(
    "102", "105", "111",
    "74"
  )) {
    # Ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    
    # Display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")
    
    # Download the data
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }
  
  
  # update the below so you remember to cite it correctly
  # knb-lter-nwt.102 (albion)
  # knb-lter-nwt.105 (gl4)
  # knb-lter-nwt.111 (martinelli)
  # knb-lter-nwt.74 (saddle)
  
  for (fname in list.files(data_dir,
                           pattern = "knb-lter.*zip",
                           full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}
# -- READ AND COMBINE SITE DATA ------------------------------------------

na_vals <- c("NaN")
data_file <- list()

data_file[["albdisch.nc.data.csv"]] <- readr::read_csv(file.path(data_dir, "albdisch.nc.data.csv"),
                                    na = na_vals, guess_max = 1000000)
data_file[["gl4disch.nc.data.csv"]] <- readr::read_csv(file.path(data_dir, "gl4disch.nc.data.csv"),
                                     na = na_vals, guess_max = 1000000)
data_file[["saddisch.nc.data.csv"]] <- readr::read_csv(file.path(data_dir, "saddisch.nc.data.csv"),
                                    na = na_vals, guess_max = 1000000)
data_file[["mardisch.nc.data.csv"]] <- readr::read_csv(file.path(data_dir, "mardisch.nc.data.csv"),
                                   na = na_vals, guess_max = 1000000)

df <- data_file %>%
  data.table::rbindlist(., use.names = TRUE, idcol = "file", fill = TRUE)

# -- INFILL MARTINELLI AND SADDLE (SEASONAL, INTERMITTENT SITES) --------
# Martinelli and Saddle both run dry for part of the year, so the
# infilling strategy differs from the perennial sites (GL4, Albion):
# known dry periods are set to zero flow directly from site notes, then
# any remaining gaps within the flow season are imputed.

# Documented dry-channel periods for Martinelli, based on site notes.
# Some of these are already reflected in the source data's own flags,
# but a few gaps could still be filled in directly from the notes here.

df <- df %>%
  # 1998 Record ends 09/24/98 when channel dry and logger removed.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "1998-09-24" & date < "1999-01-01",
  0, discharge
  )) %>%
  # 1999 Record ends when channel dry and logger removed (13:00 MST 14/09/99).
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "1999-09-14" & date < "2000-01-01",
  0, discharge
  )) %>%
  # 2012: The channel was dry after 9/08/12.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2012-08-09" & date < "2013-01-01",
  0, discharge
  )) %>%
  # 2015: The channel was dry after 10/09/15.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2015-10-15" & date < "2016-01-01",
  0, discharge
  )) %>%
  # 2015: The channel was dry after 10/09/15.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2015-10-15" & date < "2016-01-01",
  0, discharge
  )) %>%
  # 2017: The channel was dry after 21 September 2017 17:00 MST.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2017-09-17" & date < "2018-01-01",
  0, discharge
  )) %>%
  # 2018: Channel was dry after 16 August 2018 13:25 MST
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2018-08-16" & date < "2019-01-01",
  0, discharge
  )) %>%
  # 2019 Channel was dry after 21 August 2019
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2019-08-21" & date < "2020-01-01",
  0, discharge
  )) %>%
# 2020 Channel was dry after 24 Sept 2020
mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
  date >= "2020-09-24" & date < "2021-01-01",
0, discharge
))


# Determine each site's typical flow season (first/last day-of-year with
# any recorded nonzero flow), so days outside that window can be
# confidently set to zero rather than imputed.
flow_season <- df %>%
  mutate(yday = lubridate::yday(date)) %>%
  mutate(discharge = ifelse(!is.na(notes), NA, discharge)) %>%
  filter(!is.na(discharge) & discharge > 0) %>%
  group_by(local_site) %>%
  summarise(start = min(yday), end = max(yday))

# Fill the off-season (outside each site's flow window) with zero flow,
# and clear out any previously-infilled/estimated values so the
# imputation below works from raw data only.
infilled_df <- df %>%
  mutate(
    yday = lubridate::yday(date),
    year = lubridate::year(date)
  ) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
  left_join(flow_season) %>%
  # do not use any previously infilled data
  mutate(discharge = ifelse(!is.na(notes) & !notes %in% c(
    "estimate based on partial records only",
    "flow data estimated from observations",
    "flow data estimated from intermittent observation",
    "flow data estimated from intermittent observations",
    "flow data estimated from field observations"
  ), NA, discharge)) %>%
  mutate(discharge = ifelse(yday < start, 0, discharge)) %>%
  mutate(discharge = ifelse(yday > end, 0, discharge)) # %>%

# Step 2: use multivariate imputation (mnimput) to fill remaining
# within-season gaps, using the other seasonal site as a predictor plus
# all other years at the same site, restricted to the shared flow window.

# --- Saddle, imputed first ---

df_wide <- infilled_df %>%
  filter(yday > min(flow_season$start[flow_season$local_site %in% c("sdl", "mar")]) &
    yday < max(flow_season$end[flow_season$local_site %in% c("sdl", "mar")])) %>%
  filter(local_site %in% c("mar", "sdl")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
  arrange(date)

# Add each year of Saddle data as its own predictor column, so the
# imputation can draw on inter-annual similarity as well as the
# Martinelli/Saddle relationship.
for (yr in unique(df_wide$year)) {
  newdat <- df %>%
    filter(local_site == "sdl" & lubridate::year(date) == yr) %>%
    mutate(yday = lubridate::yday(date))
  if (yr == 2013) {
    newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove flood as seasonal predictor
  }
  # skip as predicators any years that are all na
  if (!all(is.na(newdat$discharge))) {
    newdat[[paste0("sdl_", yr)]] <- newdat[["discharge"]]
    df_wide <- df_wide %>%
      left_join(., newdat %>%
        select(
          -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
          -notes
        ))
  }
}

# Remove each year's self-prediction (a year cannot predict itself)
for (yr in unique(df_wide$year)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$year == yr] <- NA
  }
}

df_wide <- df_wide %>%
  # remove flood
  mutate(sdl = ifelse(year == 2013 & yday %in% c(253:259), NA, sdl)) %>%
  arrange(date)


# Multivariate EM imputation, time-series mode, spline df = sp_control

pred_sdl <- mnimput(df_wide %>%
  select(-date, -year, -yday) %>%
  as.data.frame(.),
data = df_wide %>%
  # select(-yday, -year, -alb)%>%
  select(-date, -year, -yday) %>%
  as.data.frame(.), # %>%
# tail(1000),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 3))
)


sdl_infl <- predict(pred_sdl)
sdl_infl$date <- df_wide %>% # tail(1000)%>%
  dplyr::pull(date)

sdl_infl <- sdl_infl %>%
  select(date, mar, sdl) %>%
  tidyr::pivot_longer(cols = mar:sdl, names_to = "local_site", values_to = "discharge")

# quick plot of how the infilling looks
if (plot_in_code) {
  # Visual check: do the imputed Saddle values (red) look reasonable
  # relative to observed data (blue)?
  ggplot(sdl_infl %>% filter(local_site == "sdl" &
    date >= "1998-10-01"), aes(x = date, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(color = "red") +
    geom_point(data = infilled_df %>% filter(local_site == "sdl"), aes(x = date, y = discharge), color = "blue") +
    # geom_line(color="blue")+
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}

# Bring the Saddle predictions into the main dataset. Infilling starts in
# winter 1998, since spring 1999 is the earliest year with any observed
# Saddle data.

infilled_df <- infilled_df %>%
  # only infill starting in winter of 1998, spring 1999 is the earliest year
  # with any saddle data
  left_join(sdl_infl %>% filter(local_site == "sdl" & date >= "1998-10-01") %>%
    rename(sdl_disch_inf = discharge))

if (plot_in_code) {
  # Per-water-year check of observed (blue) vs. imputed (red) values.
  # Looks reasonable overall; short (<2 day) gaps could alternatively be
  # handled with simple time-series interpolation rather than mnimput.
  ggplot(infilled_df %>% filter(local_site == "sdl" &
    date >= "1998-10-01"), aes(x = yday, y = discharge)) +
    geom_point(
      data =
        infilled_df %>% filter(local_site == "sdl" &
          date >= "1998-10-01" & is.na(discharge)),
      aes(y = sdl_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)

  ggplot(infilled_df %>% filter(local_site == "sdl"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "sdl" & is.na(discharge)),
      aes(y = sdl_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}



# --- Martinelli, imputed using the same approach, now with the
#     newly-infilled Saddle data available as a predictor ---

df_wide <- infilled_df %>%
  filter(yday > min(flow_season$start[flow_season$local_site %in% c("sdl", "mar")]) &
    yday < max(flow_season$end[flow_season$local_site %in% c("sdl", "mar")])) %>%
  filter(local_site %in% c("mar", "sdl")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
  arrange(date)

# add in the mar for all years
for (yr in unique(df_wide$year)) {
  newdat <- df %>%
    filter(local_site == "mar" & lubridate::year(date) == yr) %>%
    mutate(yday = lubridate::yday(date))
  if (yr == 2013) {
    newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
  }
  # skip as predicators any years that are all na
  if (!all(is.na(newdat$discharge))) {
    newdat[[paste0("mar_", yr)]] <- newdat[["discharge"]]
    df_wide <- df_wide %>%
      left_join(., newdat %>%
        select(
          -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
          -notes
        ))
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$year)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$year == yr] <- NA
  }
}

df_wide <- df_wide %>%
  # remove flood
  mutate(mar = ifelse(year == 2013 & yday %in% c(253:259), NA, mar)) %>%
  arrange(date)
# make predictions

pred_mar <- mnimput(df_wide %>%
  select(-date, -year, -yday) %>%
  as.data.frame(.),
data = df_wide %>%
  # select(-yday, -year, -alb)%>%
  select(-date, -year, -yday) %>%
  as.data.frame(.), # %>%
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 3))
)


mar_infl <- predict(pred_mar)
mar_infl$date <- df_wide %>% # tail(1000)%>%
  dplyr::pull(date)

mar_infl <- mar_infl %>%
  select(date, mar, sdl) %>%
  tidyr::pivot_longer(cols = mar:sdl, names_to = "local_site", values_to = "discharge")

if (plot_in_code) {
  # vis check
  ggplot(mar_infl %>% filter(local_site == "mar" & date < "2019-10-07"), aes(x = date, y = discharge)) +
    geom_line(data = infilled_df %>% filter(local_site == "mar"), aes(x = date, y = discharge), color = "red") +
    geom_point(data = infilled_df %>% filter(local_site == "mar"), aes(x = date, y = discharge), color = "red") +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}


# end infilling mart off sdl
# add to df
infilled_df <- infilled_df %>%
  left_join(mar_infl %>% filter(local_site == "mar" & date < "2019-10-07") %>%
    rename(mar_disch_inf = discharge))
if (plot_in_code) {
  # Looks reasonable overall, except 1987, where the observed record
  # appears to start too late in the season for the imputation to fill
  # in confidently -- worth a closer look in a future revision.
  ggplot(infilled_df %>% filter(local_site == "mar" & date < "2019-10-07"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "mar" &
          date < "2019-10-07" & is.na(discharge)),
      aes(y = mar_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)


  ggplot(infilled_df %>% filter(local_site == "mar" & date < "2019-10-07" & wyear == 1987), aes(x = yday, y = discharge)) +
    geom_point(
      data =
        infilled_df %>% filter(local_site == "mar" &
          date < "2019-10-07" & is.na(discharge) & wyear == 1987),
      aes(y = mar_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}



# -- INFILL GREEN LAKE 4 --------------------------------
# GL4 flows year-round, so it is infilled against its own record from
# other water years, rather than against another site -- correlation
# with other sites drops to near zero in winter months, which is
# exactly when GL4 has the most missing data.

df_wide <- infilled_df %>%
  filter(local_site %in% c("gl4")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>% # ,
  arrange(date) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1))

# add in the gl4 for all wyears
for (yr in unique(df_wide$wyear)) {
  if (yr == 2021) {
    next
  } else { # 1 day only
    if (yr == 2015) {
      next
    } else { # sensor appears bad
      newdat <- df %>%
        mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
        # month = lubridate::month(date)) %>%
        mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
        filter(local_site == "gl4" & wyear == yr)
      if (yr == 2013) {
        newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
      }

      # skip as predicators any years that are all na
      if (!all(is.na(newdat$discharge))) {
        newdat[[paste0("gl4_", yr)]] <- newdat[["discharge"]]
        df_wide <- df_wide %>%
          left_join(., newdat %>%
            select(
              -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
              -notes, -year, -wyear
            ))
      }
    }
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$wyear)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$wyear == yr] <- NA
  }
}

# Exclude wyear 2015 (bad sensor data) and the 2013 flood period from
# use as predictors, since neither is representative of typical flow.
# TODO: apply the same 2013-flood exclusion to Martinelli/Saddle above,
# for consistency.

df_wide <- df_wide %>%
  arrange(date) %>%
  mutate(gl4 = ifelse(wyear == 2015, NA, gl4)) %>%
  # sce should do this for mar, sad too
  mutate(gl4 = ifelse(wyear == 2013 & yday %in% c(253:259), NA, gl4))

# make predictions
pred_gl4 <- mnimput(df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
data = df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 4))
)

gl4_infl <- predict(pred_gl4)
gl4_infl$date <- df_wide %>%
  dplyr::pull(date)

gl4_infl <- gl4_infl %>%
  select(date, gl4) %>%
  tidyr::pivot_longer(cols = gl4, names_to = "local_site", values_to = "discharge")

if (plot_in_code) {
  # vis check
  ggplot(gl4_infl %>% filter(local_site == "gl4"), aes(x = date, y = discharge)) +
    geom_line(color = "blue") +
    geom_line(
      data = infilled_df %>% filter(local_site == "gl4"),
      aes(x = date, y = discharge), color = "red"
    ) +
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}

# join to the rest of the data
infilled_df <- infilled_df %>%
  left_join(gl4_infl %>% filter(local_site == "gl4" & date < "2020-10-01") %>%
    rename(gl4_disch_inf = discharge)) %>%
  # Exclude wyear 2015 (broken sensor) and 2018 (very little data) from
  # the final infilled values, not just the predictor set.
  mutate(
    gl4_disch_inf = ifelse(wyear %in% c(2015, 2018), NA, gl4_disch_inf),
    discharge = ifelse(wyear == 2015 & local_site == "gl4", NA, discharge)
  )

if (plot_in_code) {
  # looks okish?
  ggplot(infilled_df %>% filter(local_site == "gl4" & date < "2020-10-01"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "gl4" &
          date < "2020-10-01" & is.na(discharge)),
      aes(y = gl4_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}
# -- INFILL ALBION  --------------------------------------
# Same approach as for GL4 above.

df_wide <- infilled_df %>%
  filter(local_site %in% c("alb")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>% # ,
  arrange(date) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1))

# add in the alb for all wyears
for (yr in unique(df_wide$wyear)) {
  if (yr == 2021) {
    next # only 1 day of data available for this water year
  } else { 
    newdat <- df %>%
      mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
      mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
      filter(local_site == "alb" & wyear == yr)
    if (yr == 2013) {
      newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
    }

    # skip as predicators any years that are all na
    if (!all(is.na(newdat$discharge))) {
      newdat[[paste0("alb_", yr)]] <- newdat[["discharge"]]
      df_wide <- df_wide %>%
        left_join(., newdat %>%
          select(
            -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
            -notes, -year, -wyear
          ))
    }
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$wyear)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$wyear == yr] <- NA
  }
}

# do not use correlations with wyear 2015 in the predictions
# as these data are bad, similarly
# we expect the correlations with the flood data across years not to be relevant
df_wide <- df_wide %>%
  arrange(date) %>%
  mutate(alb = ifelse(wyear == 2013 & yday %in% c(253:259), NA, alb))

# make predictions
pred_alb <- mnimput(df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
data = df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 4))
)

alb_infl <- predict(pred_alb)
alb_infl$date <- df_wide %>%
  dplyr::pull(date)

alb_infl <- alb_infl %>%
  select(date, alb) %>%
  tidyr::pivot_longer(cols = alb, names_to = "local_site", values_to = "discharge")

# join to the rest of the data
infilled_df <- infilled_df %>%
  left_join(alb_infl %>% filter(local_site == "alb" &
    date <= max(df$date[!is.na(df$discharge) & df$local_site == "alb"])) %>%
    rename(alb_disch_inf = discharge))

if (plot_in_code) {
  # Looks reasonable overall, except 1999, where observed values are
  # consistently very low across the whole year -- possibly a
  # site-management change; worth investigating further.
  ggplot(infilled_df %>% filter(local_site == "alb" & date < "2020-10-01"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "alb" &
          date < "2020-10-01" & is.na(discharge)),
      aes(y = alb_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y")
}

# Backup copy of the data at this stage, before final assembly, in case
# any of the steps below need to be re-run or debugged separately.
bu <- infilled_df

# -- ASSEMBLE FINAL GAP-FILLED DATASET ------------------------------------

infilled_df <- infilled_df %>%
  mutate(is_infilled = ifelse(is.na(discharge), TRUE, FALSE)) %>%
  mutate(is_infilled = ifelse(yday < start | yday > end, TRUE, is_infilled)) %>%
  # Fill any remaining short gaps (<=3 days) via linear interpolation
  # before falling back on the site-specific imputed values
  group_by(local_site) %>%
  arrange(local_site, date) %>%
  mutate(discharge = zoo::na.approx(discharge, maxgap = 3, method = "linear")) %>%
  ungroup() %>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar",
    mar_disch_inf, discharge
  )) %>%
  # wyear<2020, mar_disch_inf, discharge))%>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "sdl",
    sdl_disch_inf, discharge
  )) %>% # &
  ## wyear>1998, sad_disch_inf, discharge))
  mutate(discharge = ifelse(is.na(discharge) & local_site == "gl4" & !wyear %in% c(2015, 2018),
    gl4_disch_inf, discharge
  )) %>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "alb",
    alb_disch_inf, discharge
  ))

# Re-run the short-gap linear interpolation once more, since substituting
# in the imputed values above can occasionally leave small new gaps
# (e.g. at the boundary between observed and imputed data).

infilled_df <- infilled_df %>%
  # gaps of up to 3 d in mart/sdl againn estimate by linear interpolation
  group_by(local_site) %>%
  arrange(local_site, date) %>%
  mutate(discharge = zoo::na.approx(discharge, maxgap = 3, method = "linear")) %>%
  ungroup()


# The log-scale imputation can occasionally back-transform to a handful
# of very slightly negative values; clip these to zero.

infilled_df <- infilled_df %>%
  mutate(discharge = ifelse((!is.na(discharge) & discharge < 0), 0, discharge))

# Export the gap-filled dataset for use in streamflow_data_dashboard.R
write.csv(infilled_df %>%
            select(local_site, date, discharge, yday, year, wyear, is_infilled),
          file.path(data_dir, paste0("spctl_", sp_control, ".csv")),
          row.names = FALSE
)
