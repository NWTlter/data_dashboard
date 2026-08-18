# ============================================================================
# Niwot Ridge LTER Data Dashboard: Stream Chemistry Anomaly Figures
# ============================================================================
#
# Purpose:
#   Reads September stream water sulfate chemistry data from four Niwot
#   Ridge sites and produces:
#
#     1. sulfate_conc_sept.png   September sulfate concentration (native
#                                 units, ueq/L), by site
#     2. sulfate_anom_sept.png   September sulfate anomaly, as absolute
#                                 difference from each site's long-term
#                                 September mean (native units, ueq/L --
#                                 not % change), by site
#
# Input:
#   data/albisolu.nc.data.csv   -- Albion
#   data/gre4solu.nc.data.csv   -- Green Lake 4
#   data/grrgsolu.nc.data.csv   -- Green Lake 5 Rock Glacier
#   data/ariksolu.nc.data.csv   -- Arikaree
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

# -- SETUP --------------------------------------------------------------------
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

options(stringsAsFactors = FALSE)

library(tidyverse)
library(ggthemes)  # theme_hc()
library(ggpmisc)   # symmetric_limits() for scale_y_continuous
library (EDIutils)   # read_data_package_archive(), list_data_package_revisions()
options(stringsAsFactors = F)

theme_set(theme_bw())

data_dir <- file.path("stream_water_chem","data")
figures_dir <- file.path("stream_water_chem","figures")

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

 # Authenticate with EDI using an API key stored as an environment
  # variable "EDI_API_KEY" (e.g. via .Renviron locally, or a GitHub Actions secret in
  # CI). Never hardcode the key here.
  # The key must be generated from EDI's identity and access manager
  # https://edirepository.org/resources/iam
  # it is recommended for READ use, such as this to make an access key
  # with no privedges (which will default to read), copy the key and and then set it as
  # an environment variable in your R session (Sys.setenv(EDI_API_KEY = "your_key_here"))
  # If you have already modifed your R environment to have the key set, you can skip this step.
  # If you don't have an EDI key, set download_data <- FALSE
  # and manually download the data from EDI to the /data directory.

if (download_data) {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  # Authenticate with EDI using an API key stored as an environment
  # variable "EDI_API_KEY" (e.g. via .Renviron locally, or a GitHub Actions secret in
  # CI). Never hardcode the key here.
  # The key must be generated from EDI's identity and access manager
  # https://edirepository.org/resources/iam
  # it is recommended for READ use, such as this to make an access key
  # with no privedges (which will default to read), copy the key and and then set it as
  # an environment variable in your R session (Sys.setenv(EDI_API_KEY = "your_key_here"))
  # If you have already modifed your R environment to have the key set, you can skip this step.
  # If you don't have an EDI key, set download_data <- FALSE
  # and manually download the data from EDI to the /data directory.

  scope <- "knb-lter-nwt"  # Niwot scope

  # Stream chemistry inputs used, with current citations (update as
  # package versions change):
  #
  #   Albion: 103, 1982-ongoing
  #     Caine, T., J. Morse, S. Yevak, H. Brandes, and Niwot Ridge LTER.
  #     2026. Stream water chemistry data for Albion site, 1982 - ongoing.
  #     ver 16. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/ae013a0d39b2499f3b8feda36d59b86f
  #
  #   Green Lake 4: 108, 1982-ongoing
  #     Caine, T., J. Morse, B. Gager, S. Yevak, H. Brandes, and Niwot
  #     Ridge LTER. 2026. Stream water chemistry data for Green Lake 4,
  #     1982 - ongoing. ver 15. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/4c3446c0c939769f7cef88654ecd7be4
  #
  #   Green Lake 5 Rock Glacier: 163, 1998-ongoing
  #     Caine, T., J. Morse, B. Gager, S. Yevak, H. Brandes, and Niwot
  #     Ridge LTER. 2026. Stream water chemistry data for Green Lake 5
  #     Rock Glacier, 1998 - ongoing. ver 5. Environmental Data
  #     Initiative.
  #     https://doi.org/10.6073/pasta/fb29d1c281cc934a890564f1e9bee5cf
  #
  #   Arikaree: 104, 1984-ongoing
  #     Caine, T., J. Morse, B. Gager, S. Yevak, H. Brandes, and Niwot
  #     Ridge LTER. 2026. Stream water chemistry data for Arikaree
  #     cirque, 1984 - ongoing. ver 15. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/0947b991e4e6071a0ec8f7db3e8d17db
  #

  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  for (id in c("103", "108", "163", "104")) {
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
  grrgsolu = "grrgsolu.nc.data.csv",  # Green Lake 5 Rock Glacier
  ariksolu = "ariksolu.nc.data.csv"   # Arikaree
)

data_file <- lapply(
  file.path(data_dir, site_files),
  read.csv, na = na_vals
)

datafile <- plyr::rbind.fill(data_file) %>% data.frame(.)

# -- CLEAN DATA -----------------------------------------------------------

non_numeric_cols <- c("LTER_site", "local_site", "year", "date", "time", "comments")
num_vars <- setdiff(names(datafile), non_numeric_cols)

#Treat "trace", "<mdl", and similar non-numeric flags as 0. This may not
# be the ideal treatment for every use case, but seems preferable to
# dropping those records outright.

repl.f <- function(x) ifelse(grepl("<|trace|u", x), 0, x)
datafile <- datafile %>% dplyr::mutate_at(vars(all_of(num_vars)), .funs = repl.f)


# Remove values >4 SD above that month's long-term average per site --
# there are a handful of large outliers that look like data errors on
# inspection.

flag_vals <- datafile %>%
  mutate(
    local_site = ifelse(
      local_site == "GREEN LAKE 4 WATERFALL", "GREEN LAKE 4", local_site
    ),
    `SO4..` = as.numeric(as.character(`SO4..`)),
    month = lubridate::month(date)
  ) %>%
  group_by(local_site, month) %>%
  summarise(
    sulfate = mean(SO4.., na.rm = TRUE),
    sulfate_sd = sd(SO4.., na.rm = TRUE),
    sulfate_cutoff = sulfate + 4 * sulfate_sd,
    .groups = "drop"
  ) %>%
  select(local_site, month, sulfate_cutoff)

# The four sites this script focuses on -- Albion, Green Lake 4, Green
# Lake 5 Rock Glacier, and Arikaree. Saddle Stream, Martinelli, and Green
# Lake 5 Outlet are intentionally excluded because of depracated sampling
# GL5 outlet) and/or no flows in sept (saddle, martinelli).
sites_used <- c("ALBION", "GREEN LAKE 4", "GREEN LAKE 5 ROCK GLACIER", "ARIKAREE")

# ============================================================================
# FIGURE 1: September sulfate concentration (native units), by site
# ============================================================================
# Late-season (September) sulfate is used here as an indicator of
# late-summer sulfur flux; concentrations are plotted in their native
# units (ueq/L), with each site's long-term September mean shown as a
# reference line.

chem_by_month <- datafile %>%
  mutate(
    local_site = ifelse(
      local_site == "GREEN LAKE 4 WATERFALL", "GREEN LAKE 4", local_site
    ),
    `SO4..` = as.numeric(as.character(`SO4..`)),
    month = lubridate::month(date)
  ) %>%
  left_join(., flag_vals) %>%
  mutate(
    `SO4..` = ifelse(`SO4..` > sulfate_cutoff, NA, `SO4..`)
  ) %>%
  group_by(year, local_site, month) %>%
  summarise(
    sulfate = mean(SO4.., na.rm = TRUE),
    .groups = "drop"
  )


sulfate_sept <- chem_by_month %>%
  filter(month == 9 & local_site %in% sites_used) %>%
  group_by(local_site) %>%
  mutate(
    mean_S = mean(sulfate, na.rm = TRUE)
  ) %>%
  mutate(
    anom_S_abs = sulfate - mean_S,
    posneg = ifelse(anom_S_abs > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  ) %>%
  ungroup()

x_range_S <- range(sulfate_sept$year, na.rm = TRUE)

sulfate_y_label <- expression(SO[4]^{"2-"} ~ "(" * mu * "Eq " * L^-1 * "; September)")

g1 <- ggplot(sulfate_sept, aes(x = year, y = sulfate)) +
  geom_col(fill = "#8C2F12") +
  labs(y = sulfate_y_label, x = "Year") +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(sulfate_sept$year) / 10),
      10 * ceiling(max(sulfate_sept$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(sulfate_sept$year) / 5),
      5 * ceiling(max(sulfate_sept$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  theme_hc() +
  facet_wrap(~local_site, scales = "free_y") +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5)
  )

ggsave(g1,
       file = file.path(figures_dir, "sulfate_conc_sept.png"),
       scale = 0.5, width = 15, height = 6
)

# ============================================================================
# FIGURE 2: September sulfate anomaly (absolute, native units), by site
# ============================================================================
# Same underlying data as Figure 1, but expressed as the absolute
# difference from each site's long-term September mean (ueq/L), rather
# than a percent change.

row_labels_S <- data.frame(
  # facet_wrap sorts alphabetically: row 1 = Albion, Arikaree; row 2 =
  # Green Lake 4, Green Lake 5 Rock Glacier -- these are each row's
  # rightmost facet.
  local_site = c("ARIKAREE", "GREEN LAKE 5 ROCK GLACIER"),
  x = x_range_S[2] + diff(x_range_S) * 0.1,
  y = 0,  # each facet's y-scale is symmetric around 0, so center the label vertically
  label = "more~sulfate  %<->%  less~sulfate"
)

# scale_y_continuous(limits = symmetric_limits) only computes one global
# symmetric range -- with free_y facets we instead need each site
# individually centered on 0, so supply per-site +/- limits via an
# invisible geom_blank() layer (one row per site, at its own max
# absolute anomaly).
symmetric_limits_by_site <- sulfate_sept %>%
  group_by(local_site) %>%
  summarise(max_abs = max(abs(anom_S_abs), na.rm = TRUE), .groups = "drop") %>%
  tidyr::crossing(sign = c(-1, 1)) %>%
  mutate(anom_S_abs = max_abs * sign)

g2 <- ggplot(sulfate_sept, aes(x = year, y = anom_S_abs)) +
  geom_col(aes(fill = posneg)) +
  geom_blank(data = symmetric_limits_by_site, aes(y = anom_S_abs), inherit.aes = FALSE) +
  scale_fill_manual(values = c("#8C2F12", "#F5DDB0")) +
  labs(
    y = expression(atop("September sulfate anomaly (" * mu * "Eq " * L^-1 * ")", "difference from long-term mean")),
    x = "Year"
  ) +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(sulfate_sept$year) / 10),
      10 * ceiling(max(sulfate_sept$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(sulfate_sept$year) / 5),
      5 * ceiling(max(sulfate_sept$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  geom_text(
    data = row_labels_S,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  facet_wrap(~local_site, scales = "free_y") +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)  # extra right margin for the row label
  ) +
  coord_cartesian(xlim = x_range_S, clip = "off")

ggsave(g2,
       file = file.path(figures_dir, "sulfate_anom_sept.png"),
       scale = 0.5, width = 15, height = 6
)
