# ============================================================================
# Niwot Ridge LTER Data Dashboard: Pika Recruitment Figure
# ============================================================================
#
# Purpose:
#   Reads pika capture data (juveniles vs. adults) and produces:
#
#     1. pika_ja_ratio.png   Annual juvenile:adult capture ratio anomaly,
#                             expressed as a log2 ratio relative to the
#                             long-term mean, at West Knoll
#
# Input:
#   data/pikas_no_recap.cr.data.csv     -- historical captures, 1981-1990
#                                           (EDI knb-lter-nwt.43)
#   data/pika_demography.cr.data.csv    -- current captures, 2008-ongoing
#                                           (EDI knb-lter-nwt.8, West Knoll only)
#
# Output:
#   One PNG figure written to plot_dir (see file name above).
#
# Original analysis: Sarah C. Elmendorf, 12 Jan 2022 (based on an earlier
#   script by Chris Ray: script-pika-juvie-capture-ratio-vs-GDD.R)
# Updated to use log ratios: Sarah C. Elmendorf, 6 July 2026
# Revised and updated for the Data Dashboard: Anne Marie Panetta, July 2026
# ============================================================================

# -- SETUP --------------------------------------------------------------
# libraries
library(tidyverse)
library(EDIutils) # tools for interacting with EDI's data package API
library (scales)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


data_dir <- file.path("data")
plot_dir <- file.path("plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Set to TRUE to (re)download source data from EDI. Only needs to be run
# once, or when a newer package version is released.
download_data <- FALSE

# download data -----------------------------------------------------------
# Note: if you have already downloaded SOME data, the read_data_package_archive()
#  will overwrite existing files. Clear your /data directory first if you need
# to forse a fresh download.

if (download_data) {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  scope <- "knb-lter-nwt" # Niwot scope
  
  # Data packages:
  #   43 -- small mammal species composition data for Niwot Ridge, 1981-1990
  #         (includes the pika-cleared-of-recapture table)
  #     Halfpenny, J., C. Ray, and Niwot Ridge LTER. 2025. Small mammal
  #     species composition data for Niwot Ridge, 1981 - 1990. ver 6.
  #     Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/129b75c156d4e825e691b6d5ed401300
  #
  #   8  -- pika demography data for west knoll and Indian Peaks wilderness,
  #         2008-ongoing
  #     Ray, C. and Niwot Ridge LTER. 2025. Pika demography data for west
  #     knoll and Indian Peaks wilderness, 2008 - ongoing. ver 9.
  #     Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/702d00888da15e13ce9a6daf36dcc06d
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  
  for (id in c("43", "8")) {
    # Ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")

    # Display current version -> this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")

    # Download the data
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }

  # Unzip all zip files in the directory
  # This overwrites the manifests, but you don't really need them.
  for (fname in list.files(data_dir,
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}

# -- READ DATA ------
# "d" : file (pika data cleared of recaptures) from EDI package knb-lter-nwt.43
# Note: this table only spans 1981-1990 (year-unique captures), unlike the
# previously-used local pika-demog-WK-1981-2020.csv which combined this table
# with more recent demography data

d <- read_csv(file.path(data_dir, "pikas_no_recap.cr.data.csv"))

# prep pika data
d <- mutate(d,
  date = as.Date(date),
  day = lubridate::mday(date),
  month = lubridate::month(date),
  year = lubridate::year(date),
  doy = lubridate::yday(date),
  count = 1
)

# "d2" file (pika demography, 2008-ongoing) from EDI package knb-lter-nwt.8
# Restricted to west knoll (WK) for continuity with the historical "d" series
# (pkg 8 also covers WKN, LL, ML, CG, WB, GLV4, D1, added at various points)

d2 <- read_csv(file.path(data_dir, "pika_demography.cr.data.csv")) %>%
  mutate(weight = ifelse(weight == '159?', "159", weight)) %>%
  mutate(weight = as.numeric(weight)) %>%
  filter(local_site == "WK") %>%
  mutate(
    date = as.Date(date),
    day = lubridate::mday(date),
    month = lubridate::month(date),
    year = lubridate::year(date),
    doy = lubridate::yday(date),
    count = 1,
    # build a tag ID from both ears (either can be lost/damaged), matching
    # the id_num convention already used in "d"
    id_num = paste0(code_r_ear, num_r_ear, code_l_ear, num_l_ear)
  ) %>%
  group_by(id_num, year) %>%
  mutate(
    # if this capture's stage is uncertain/unsampled (A?, J?, NS, NA), borrow
    # the confirmed stage from another capture of the same individual in the
    # same year, if one exists; drop the individual-year if never confirmed
    stage = case_when(
      any(stage == "A") ~ "A",
      any(stage == "J") ~ "J",
      any(stage == "A?") ~ "A",
      any(stage == "J?") ~ "J",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(stage)) %>%
  # one record per individual per year: keep the first capture of the year
  arrange(id_num, year, date) %>%
  distinct(id_num, year, .keep_all = TRUE) %>%
  select(date, id_num, sex, stage, weight, year, month, day, doy, count)

# combine the historical (pkg 43) and current (pkg 8) pika tables
d_all <- bind_rows(
  d %>% select(date, id_num, sex, stage, weight, year, month, day, doy, count),
  d2
)

by.yr <- summarize(group_by(d_all, year),
  juvies = sum(count[stage == "J"]),
  adults = sum(count[stage == "A"]),
  capdate.mu = mean(doy, na.rm = T)
)


# -- DERIVE ANOMALIES -------------------------------------------------------
by.yr.anom <- by.yr %>%
  mutate(j.a = juvies / adults) %>%
  mutate(
    j.a = ifelse(!is.finite(j.a), NA, j.a),
    # log-ratio scale: anomalies are relative to the mean log-ratio, so
    # juvenile-heavy and adult-heavy years are symmetric around zero
    log_ratio = log2(j.a),
    log_ratio = ifelse(!is.finite(log_ratio), NA, log_ratio),
    mean_log_ratio = mean(log_ratio, na.rm = TRUE),
    posneg_log = ifelse(log_ratio > mean_log_ratio, "pos", "neg") %>%
      factor(c("pos", "neg")),
  )

# ============================================================================
# FIGURE 1: Annual pika juvenile:adult capture ratio anomalies at West Knoll
# ============================================================================
# Everything is plotted on a log-ratio scale: bars originate at the mean
# log-ratio (dashed grey line) rather than at 0, so years above the mean
# are purple and rise, years below are gold and fall. Axis ticks/labels
# are back-transformed from log-ratio to ratio units. A dotted black line
# marks the 1:1 point (equal juveniles and adults captured) as a fixed,
# biologically meaningful reference distinct from the dataset's own mean.

bar_half_width <- 0.4

# Symmetric y-limits around the 1:1 line (log2(1) = 0), so juvenile-heavy and
# adult-heavy years get equal visual space
y_max_abs <- max(abs(range(by.yr.anom$log_ratio, na.rm = TRUE)))

# Fixed x-range for the panel itself (bar extents only) - kept explicit so
# the panel does NOT auto-expand to include the label annotation below
year_range <- diff(range(by.yr.anom$year, na.rm = TRUE))
x_panel_lim <- range(by.yr.anom$year, na.rm = TRUE) + c(-bar_half_width, bar_half_width)

# x position for the manual secondary-axis label: comfortably past the
# fixed panel edge (and past the right-hand tick marks), in the blank margin
x_label_pos <- max(by.yr.anom$year, na.rm = TRUE) + year_range * 0.15

# Match the annotation's font to the left axis title's actual size (theme_hc's
# axis.title.y, in pt); annotate()/geom_text size is in mm, hence the /.pt
axis_title_size <- calc_element("axis.title.y", theme_hc())$size / .pt

# Secondary-axis label placed manually (rather than via sec_axis(name = ...))
# so it can be centered on the 1:1 (log2(1) = 0) break -- ggplot always
# vertically centers a sec_axis name on the whole panel instead of on a
# specific break.
x_label_pos <- max(by.yr.anom$year, na.rm = TRUE) + year_range * 0.15

row_label <- data.frame(
  x = x_label_pos,
  y = 0,
  label = "higher~recruitment  %<->%  lower~recruitment"  # plotmath expression for the arrow glyph
)

g1 <- ggplot(by.yr.anom) +
  geom_rect(aes(
    xmin = year - bar_half_width, xmax = year + bar_half_width,
    ymin = pmin(mean_log_ratio, log_ratio), ymax = pmax(mean_log_ratio, log_ratio),
    fill = posneg_log
  )) +
  scale_fill_manual(values = c("#5E3C99", "#E1A100")) +
  geom_hline(
    yintercept = unique(by.yr.anom$mean_log_ratio),
    linetype = "longdash", color = "gray70", linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted", color = "black", linewidth = 0.6
  ) +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(by.yr.anom$year) / 10),
      10 * ceiling(max(by.yr.anom$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(by.yr.anom$year) / 5),
      5 * ceiling(max(by.yr.anom$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)  # 5-yr minor ticks alongside 10-yr labeled major ticks
  ) +
  scale_y_continuous(
    breaks = log2(c(1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8)),
    labels = c("1:8", "1:4", "1:2", "1:1", "2:1", "4:1", "8:1"),
    sec.axis = sec_axis(trans = I, breaks = 0, labels = "")
  ) +
  labs(y = "Pika juvenile:adult ratio", x = "Year") +
  theme_hc() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "gray50"),
    axis.minor.ticks.length = rel(0.5),
    plot.title = element_text(margin = margin(b = 0)),
    plot.margin = margin(t = 5.5, r = 44, b = 5.5, l = 5.5)
  ) +
  # secondary-axis title placed manually so it's centered on the 1:1 (log2(1))
  # break, since ggplot always vertically centers a sec_axis name on the
  # whole panel instead of on a specific break
  geom_text(
    data = row_label,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = axis_title_size,
    inherit.aes = FALSE, parse = TRUE
  ) +
  coord_cartesian(
    xlim = x_panel_lim,
    ylim = c(-y_max_abs, y_max_abs),
    clip = "off"
  )

ggsave(g1,
  file = file.path(plot_dir, "pika_ja_ratio.png"),
  scale = 0.6, width = 10, height = 6
)
