# ============================================================================
# Niwot Ridge LTER Data Dashboard: Stream Chemistry Anomaly Figures
# ============================================================================
#
# Purpose:
#   Reads stream water chemistry data from six Niwot Ridge sites and
#   produces:
#
#     1. sulfate_anom_sept.png     September sulfate anomaly, by site
#     2. nitrate_anom_summer.png   Summer (May-Sep) nitrate anomaly, by site
#
# Input:
#   data/albisolu.nc.data.csv   -- Albion
#   data/gre4solu.nc.data.csv   -- Green Lake 4
#   data/saddsolu.nc.data.csv   -- Saddle Stream
#   data/martsolu.nc.data.csv   -- Martinelli basin
#   data/grrgsolu.nc.data.csv   -- Green Lake 5 Rock Glacier
#   data/gre5solu.nc.data.csv   -- Green Lake 5 Outlet
#
# Output:
#   Two PNG figures written to figures_dir (see file names above).
#
# Notes:
#   See the nwt-8 long-term trends repo for a flow-normalized workflow, if
#   that's more appropriate for your question than the raw concentrations
#   used here.
#
# Original analysis: Sarah C. Elmendorf, 2023
# Revised and updated for the Data Dashboard: Anne Marie Panetta, 2026
# ============================================================================

# -- SETUP --------------------------------------------------------------

rm(list = ls())


library(tidyverse)
library(ggthemes)  # theme_hc()
library(ggpmisc)   # symmetric_limits() for scale_y_continuous
options(stringsAsFactors = F)

theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_dir <- "data"
figures_dir <- "figures"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Set to TRUE to (re)download source data from EDI. Only needs to be run
# once, or when a newer package version is released.
download_data <- TRUE

# -- DOWNLOAD DATA FROM EDI ---------------------------------------------
# Note: if you have already downloaded SOME data, read_data_package_archive()
# will not overwrite existing files. Clear the /data directory first if you
# need to force a fresh download.

if (download_data) {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  # Authenticate with EDI using an API key stored as an environment
  # variable (e.g. via .Renviron locally, or a GitHub Actions secret in
  # CI). Never hardcode the key here.
  EDIutils::login(key = Sys.getenv("EDI_API_KEY"))

  scope <- "knb-lter-nwt"  # Niwot scope
  
  # Discharge/chemistry inputs used, with current citations (update as
  # package versions change):
  #
  #   Albion: 103, 1982-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2025. Stream water chemistry data for Albion site, 1982 - ongoing.
  #     ver 15. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/39b37d16657733f1d661de4c727568a3
  #
  #   Green Lake 4: 108, 1982-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2025. Stream water chemistry data for Green Lake 4, 1982 - ongoing.
  #     ver 14. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/cdbd25e2ce7a8b17f785bd4f07aba4fe
  #
  #   Martinelli basin: 112, 1984-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2025. Stream water chemistry data for Martinelli basin,
  #     1984 - ongoing. ver 5. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/4a761ed1f7b24d25c3eb870662213789
  #
  #   Green Lake 5 Rock Glacier: 163, 1998-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2025. Stream water chemistry data for Green Lake 5 Rock Glacier,
  #     1998 - ongoing. ver 4. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/990a451879638d17eda4526215114390
  #
  #   Green Lake 5 Outlet: 109, 1984-ongoing (through ~2021)
  #     Caine, T. 2021. Stream water chemistry data for Green Lake 5
  #     outlet, 1984 - ongoing. ver 12. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/effe6d9f1a764ebaccb7ae9b28cba7a5
  #
  #   Saddle Stream site: 160, 1994-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2025. Stream water chemistry data for Saddle Stream site,
  #     1994 - ongoing. ver 4. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/92e6f11e89b730321df0865bec4b7b0b
  #
  # Other packages available in this series if needed: 162 (watershed
  # flume), 213 (Soddie).
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  for (id in c("103", "108", "160", "112", "163", "109")) {
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    packageID <- paste(scope, id, revision, sep = ".")
    
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID)) # confirm you're citing the version actually used
  }
  
  for (fname in list.files(data_dir, pattern = "knb-lter.*zip", full.names = TRUE)) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}


# -- READ DATA ------------------------------------------------------------

#  NA codes used across these files: not-collected, not-sampled, below
# detection, etc. -- all treated as missing-completely-at-random (MCAR).

na_vals <- c("", "NA", "NaN", "NP", "DNS", "NSS", "EQCL", "QNS", "NV", "dns")

site_files <- c(
  albisolu = "albisolu.nc.data.csv",  # Albion
  gre4solu = "gre4solu.nc.data.csv",  # Green Lake 4
  saddsolu = "saddsolu.nc.data.csv",  # Saddle Stream
  martsolu = "martsolu.nc.data.csv",  # Martinelli basin
  grrgsolu = "grrgsolu.nc.data.csv",  # Green Lake 5 Rock Glacier
  gre5solu = "gre5solu.nc.data.csv"   # Green Lake 5 Outlet
)

data_file <- lapply(
  file.path(data_dir, site_files),
  read.csv, na = na_vals
)

datafile <- plyr::rbind.fill(data_file) %>% data.frame(.)

# -- CLEAN DATA -----------------------------------------------------------

num_vars <- names(datafile)[6:44]

#Treat "trace", "<mdl", and similar non-numeric flags as 0. This may not
# be the ideal treatment for every use case, but seems preferable to
# dropping those records outright.

repl.f <- function(x) ifelse(grepl("<|trace|u", x), 0, x)
datafile <- datafile %>% dplyr::mutate_at(vars(num_vars), .funs = repl.f)


# Remove values >3 SD above that month's long-term average per site --
# there are a handful of large outliers that look like data errors on
# inspection.

flag_vals <- datafile %>%
  mutate(
    local_site = ifelse(
      local_site == "GREEN LAKE 4 WATERFALL", "GREEN LAKE 4", local_site
    ),
    `NO3.` = as.numeric(as.character(`NO3.`)),
    `SO4..` = as.numeric(as.character(`SO4..`)),
    month = lubridate::month(date)
  ) %>%
  group_by(local_site, month) %>%
  summarise(
    nitrate = mean(NO3., na.rm = TRUE),
    sulfate = mean(SO4.., na.rm = TRUE),
    nitrate_sd = sd(NO3., na.rm = TRUE),
    sulfate_sd = sd(SO4.., na.rm = TRUE),
    n_cutoff = nitrate + 3 * nitrate_sd,
    sulfate_cutoff = sulfate + 3 * sulfate_sd,
    .groups = "drop"
  ) %>%
  select(local_site, month, n_cutoff, sulfate_cutoff)

# ============================================================================
# FIGURE 1: September sulfate anomaly, by site
# ============================================================================
# Late-season (September) sulfate is used here as an indicator of
# late-summer sulfur flux; anomalies are percent difference from each
# site's long-term September mean.

chem_by_month <- datafile %>%
  mutate(
    local_site = ifelse(
      local_site == "GREEN LAKE 4 WATERFALL", "GREEN LAKE 4", local_site
    ),
    `NO3.` = as.numeric(as.character(`NO3.`)),
    `SO4..` = as.numeric(as.character(`SO4..`)),
    month = lubridate::month(date)
  ) %>%
  left_join(., flag_vals) %>%
  mutate(
    `NO3.` = ifelse(`NO3.` > n_cutoff, NA, `NO3.`),
    `SO4..` = ifelse(`SO4..` > sulfate_cutoff, NA, `SO4..`)
  ) %>%
  group_by(year, local_site, month) %>%
  summarise(
    nitrate = mean(NO3., na.rm = TRUE),
    sulfate = mean(SO4.., na.rm = TRUE),
    .groups = "drop"
  )


anom_chem_month <- chem_by_month %>%
  filter(month == 9 & local_site %in% c(
    "ALBION", "GREEN LAKE 4", "GREEN LAKE 5",
    "MARTINELLI"
  )) %>%
  group_by(local_site) %>%
  mutate(
    mean_S = mean(sulfate, .groups = "drop")
  ) %>%
  mutate(
    anom_S = sulfate * 100 / mean_S - 100,
    posneg = ifelse(anom_S > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  ) %>%
  ungroup()

row_labels_S <- data.frame(
  local_site = c("GREEN LAKE 4", "MARTINELLI"),
  x = x_range_S[2] + diff(x_range_S) * 0.1,
  y = 0.03 * diff(range(anom_chem_month$anom_S, na.rm = TRUE)),  # small upward nudge to visually center on 0
  label = "more~sulfate  %<->%  less~sulfate"
)

g1 <- ggplot(anom_chem_month, aes(x = year, y = anom_S)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#8C2F12", "#F5DDB0")) +
  labs(y = "September sulfate anomaly \n (% difference from long-term mean)", x = "Year") +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(anom_chem_month$year) / 10),
      10 * ceiling(max(anom_chem_month$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(anom_chem_month$year) / 5),
      5 * ceiling(max(anom_chem_month$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_continuous(limits = symmetric_limits) +
  geom_text(
    data = row_labels_S,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  facet_wrap(~local_site) +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)  # extra right margin for the row label
  ) +
  coord_cartesian(xlim = x_range_S, clip = "off")

ggsave(g1,
       file = "figures/sulfate_anom_sept.png",
       scale = 0.5, width = 15, height = 6
)

# ============================================================================
# FIGURE 2: Summer (May-Sep) nitrate anomaly, by site
# ============================================================================
# Growing-season nitrate, averaged May through September; anomalies are
# percent difference from each site's long-term summer mean.

chem_by_summer <- datafile %>%
  mutate(
    local_site = ifelse(
      local_site == "GREEN LAKE 4 WATERFALL", "GREEN LAKE 4", local_site
    ),
    `NO3.` = as.numeric(as.character(`NO3.`)),
    `SO4..` = as.numeric(as.character(`SO4..`)),
    month = lubridate::month(date)
  ) %>%
  left_join(., flag_vals) %>%
  mutate(
    `NO3.` = ifelse(`NO3.` > n_cutoff, NA, `NO3.`),
    `SO4..` = ifelse(`SO4..` > sulfate_cutoff, NA, `SO4..`)
  ) %>%
  filter(month %in% c(5:9)) %>%
  group_by(year, local_site) %>%
  summarise(
    nitrate = mean(NO3., na.rm = TRUE),
    sulfate = mean(SO4.., na.rm = TRUE),
    .groups = "drop"
  )

anom_chem_summer_N <- chem_by_summer %>%
  filter(local_site %in% c("ALBION", "GREEN LAKE 4", "GREEN LAKE 5", "MARTINELLI")) %>%
  group_by(local_site) %>%
  mutate(
    avg_N = mean(nitrate, na.rm = TRUE, .groups = "drop")
  ) %>%
  mutate(
    anom_N = nitrate * 100 / avg_N - 100,
    posneg = ifelse(anom_N > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  ) %>%
  ungroup()

x_range_N <- range(anom_chem_summer_N$year, na.rm = TRUE)
row_labels_N <- data.frame(
  local_site = c("GREEN LAKE 4", "MARTINELLI"),  # rightmost facet of each row
  x = x_range_N[2] + diff(x_range_N) * 0.1,
  y = 0.03 * diff(range(anom_chem_summer_N$anom_N, na.rm = TRUE)),  # small upward nudge to visually center on 0
  label = "more~nitrate  %<->%  less~nitrate"
)

g2 <- ggplot(anom_chem_summer_N, aes(x = year, y = anom_N)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#1B7A6F", "#8ED1C4")) +
  labs(y = "Summer nitrate anomaly \n (% difference from long-term mean)", x = "Year") +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(anom_chem_summer_N$year) / 10),
      10 * ceiling(max(anom_chem_summer_N$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(anom_chem_summer_N$year) / 5),
      5 * ceiling(max(anom_chem_summer_N$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_continuous(limits = symmetric_limits) +
  geom_text(
    data = row_labels_N,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  facet_wrap(~local_site) +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)
  ) +
  coord_cartesian(xlim = x_range_N, clip = "off")

ggsave(g2,
       file = "figures/nitrate_anom_summer.png",
       scale = 0.5, width = 15, height = 6
)
