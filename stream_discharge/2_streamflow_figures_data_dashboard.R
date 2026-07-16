# ============================================================================
# Niwot Ridge LTER Data Dashboard: Streamflow Figures
# ============================================================================
#
# Purpose:
#   Reads gap-filled daily discharge data (produced by 1_streamflow_infill_and_q_50_data_dashboard.R)
#   for four Niwot Ridge fiel sites-- Albion (alb), Green Lake 4 (gl4),
#   Martinelli (mar), and Saddle (sdl) -- and produces three summary figures:
#
#     1. streamflow_anom.png            Annual streamflow anomaly by site
#     2. streamflow_anom_gl4_month.png  Monthly (May-Sep) streamflow anomaly
#                                        at Green Lake 4
#     3. streamflow_summer_gl4_2024.png Daily 2024 summer streamflow at
#                                        Green Lake 4 vs. historical range
#
# Input:
#   data/spctl_<sp_control>.csv  -- gap-filled discharge data, one row per
#                                    site-day, with columns:
#                                    local_site, date, discharge, yday, year,
#                                    wyear, is_infilled
#
# Output:
#   Three PNG figures written to figures_dir (see file names above).
#
# Original analysis: Sarah C. Elmendorf, 10 Jan 2023
# Revised and updated for the Data Dashboard: Anne Marie Panetta, 2026 
# ============================================================================


# -- SETUP --------------------------------------------------------------

rm(list = ls())

library(tidyverse)   # data wrangling + ggplot2
library(lemon)       # facet/legend helpers
library(ggthemes)    # theme_hc()
library(lubridate)   # date handling
library(ggplotFL)    # geom_flquantiles(), for the historical percentile band
library(ggtext)      # element_markdown(), for styled axis titles

# Anchor working directory to this script's location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_dir <- "data"
figures_dir <- "figures"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# The "spline control" (sp_control) parameter used during gap-filling;
# see 1_streamflow_infill_and_q_50_data_dashboard.R for details. It has little effect on the final
# infilled values, but is included in the filename to keep multiple
# infilling runs distinguishable.
sp_control <- 7

infilled_df <- read_csv(file.path(data_dir, paste0("spctl_", sp_control, ".csv")))


# -- DERIVED METRICS ------------------------------------------------------
# The following blocks compute summary streamflow statistics (proportion
# of each year that was infilled, timing of peak flow, and flow-timing
# quantiles q20/q50/q80) that support data-quality review and exploratory
# analysis. They are not required for the three saved figures below, but
# are kept here since later dashboard iterations may draw on them.

# Proportion of each site-year that was infilled (data-quality flag),
# weighted by each day's typical share of annual discharge
wts <- infilled_df %>%
  filter(is_infilled == FALSE) %>%
  group_by(local_site, yday) %>%
  summarize(mean_daily = mean(discharge, na.rm = TRUE)) %>%
  left_join(
    infilled_df %>%
      filter(is_infilled == FALSE) %>%
      group_by(local_site, yday) %>%
      summarize(mean_daily = mean(discharge, na.rm = TRUE)) %>%
      filter(yday != 366) %>%
      group_by(local_site) %>%
      summarize(tot = sum(mean_daily))
  ) %>%
  mutate(prop = mean_daily / tot) %>%
  full_join(., infilled_df) %>%
  filter(yday < 366) %>%  # drop leap day so all years have equal-length totals
  group_by(local_site, wyear, is_infilled) %>%
  summarise(tot_disch = sum(prop, na.rm = TRUE)) %>%
  pivot_wider(names_from = is_infilled, values_from = tot_disch) %>%
  mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(prop_infilled = `TRUE` / (`FALSE` + `TRUE`))

# Timing of peak flow each water year (2-day rolling mean, to reduce noise
# from single-day sensor spikes)
peak_flow <- infilled_df %>%
  group_by(wyear, local_site) %>%
  mutate(disch_ma_2 = slider::slide_dbl(discharge, mean, .before = 2, .after = 2)) %>%
  summarize(max_disch = max(disch_ma_2, na.rm = TRUE)) %>%
  full_join(infilled_df) %>%
  mutate(disch_ma_2 = slider::slide_dbl(discharge, mean, .before = 2, .after = 2)) %>%
  rowwise() %>%
  filter(max_disch == disch_ma_2) %>%
  filter(wyear < 2025)

# Note: Albion winter floods (and the 2013 flood more broadly) can distort
# peak-timing patterns and may warrant exclusion in future refinements.
# Overall, peak timing shows less trend than the flow-timing quantiles below
# (e.g. q50), which do trend earlier over the record.
ggplot(peak_flow, aes(x = wyear, y = yday)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site) +
  ggtitle("Timing of peak flow by water year")

# Annual total discharge and flow-timing thresholds (20th/50th/80th
# percentile of cumulative annual discharge), calculated per water year
# and, separately, against the site's all-time average (for cross-year
# comparability)
tots <- infilled_df %>%
  filter(yday < 366) %>%
  group_by(local_site, wyear) %>%
  summarise(
    tot_disch = sum(discharge, na.rm = TRUE),
    tot_day = sum(!is.na(discharge))
  ) %>%
  mutate(tot_disch = ifelse(tot_day != 365, NA, tot_disch)) %>%  # only use complete years
  mutate(
    disch_50 = tot_disch / 2,
    disch_20 = tot_disch / 5,
    disch_80 = 4 * (tot_disch / 5)
  ) %>%
  left_join(., wts %>% select(local_site, wyear, prop_infilled))

avg_tots <- tots %>%
  group_by(local_site) %>%
  summarize(tot_disch = mean(tot_disch, na.rm = TRUE)) %>%
  mutate(
    disch_50_all = tot_disch / 2,
    disch_20_all = tot_disch / 5,
    disch_80_all = 4 * (tot_disch / 5)
  )

# Helper: day-of-year on which cumulative discharge first exceeds a given
# threshold. Repeated below for the 20th/50th/80th percentile thresholds,
# both per-water-year (q20/q50/q80) and against the all-time average
# (q20_mean/q50_mean/q80_mean).
q50 <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_50, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q50 = min(yday), .groups = "drop") %>%
  left_join(., tots)

q50_mean <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(avg_tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_50_all, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q50_mean = min(yday), .groups = "drop")

q20 <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_20, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q20 = min(yday), .groups = "drop")

# Note: q20 is currently always positive (day-of-year within the water
# year); for sites like Albion where 20% cumulative discharge could
# plausibly occur in the prior fall, this may need to allow negative
# values in a future revision.
q20_mean <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(avg_tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_20_all, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q20_mean = min(yday), .groups = "drop")

q80 <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_80, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q80 = min(yday), .groups = "drop")

q80_mean <- infilled_df %>%
  group_by(wyear, local_site) %>%
  arrange(date) %>%
  mutate(cum_discharge = cumsum(discharge)) %>%
  ungroup() %>%
  left_join(avg_tots) %>%
  rowwise() %>%
  mutate(thresh = ifelse(cum_discharge > disch_80_all, 1, 0)) %>%
  group_by(wyear, local_site) %>%
  filter(thresh == 1) %>%
  summarize(q80_mean = min(yday), .groups = "drop")

q50 <- q50 %>%
  left_join(., q20) %>%
  left_join(., q80) %>%
  left_join(., q50_mean) %>%
  left_join(., q20_mean) %>%
  left_join(., q80_mean) %>%
  mutate(
    duration = q80 - q20,            # spread of the middle 60% of flow, per water year
    duration_mean = q80_mean - q20_mean  # same, relative to all-time average
  )

ggplot(q50, aes(x = wyear, y = q50, color = prop_infilled)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site) +
  ggtitle("Day-of-year reaching 50% of annual discharge, over time")

# Data completeness check: confirms each site-water-year has the expected
# number of daily records. All years are complete except Albion 1981
# (partial record start) and a handful of Saddle days in the final year
# after flow had already reached zero.
cts <- infilled_df %>%
  group_by(local_site, wyear) %>%
  summarize(n = dplyr::n())


# ============================================================================
# FIGURE 1: Annual streamflow anomaly by site
# ============================================================================
# Percent difference between each water year's total discharge and that
# site's long-term average, colored by direction (above/below average).

anom_disch <- infilled_df %>%
  filter(wyear > 1981 & wyear < 2025) %>%
  group_by(local_site, wyear) %>%
  summarise(
    total_disch = sum(discharge),
    ct = dplyr::n()
  ) %>%
  ungroup() %>%
  group_by(local_site) %>%
  mutate(
    avgdisch = mean(total_disch, na.rm = TRUE),
    anom_disch = total_disch * 100 / avgdisch - 100,
    posneg = ifelse(anom_disch > 0, "pos", "neg") %>% factor(c("pos", "neg"))
  ) %>%
  ungroup()

facet_names <- c(
  alb = "Albion",
  gl4 = "Green Lake 4",
  mar = "Martinelli",
  sdl = "Saddle"
)

# "more flow / less flow" annotation is drawn once per facet row (on the
# rightmost panel), rather than once per panel, to avoid visual clutter.
# In this 2x2 grid that's Green Lake 4 (row 1) and Saddle (row 2).
row_label_sites <- c("gl4", "sdl")

x_range <- range(anom_disch$wyear, na.rm = TRUE)

row_labels <- data.frame(
  local_site = row_label_sites,
  x = x_range[2] + diff(x_range) * 0.1,  # position just past the last data year
  y = 0,
  label = "more~flow  %<->%  less~flow"  # plotmath expression for the arrow glyph
)

g1 <- ggplot(anom_disch, aes(x = wyear, y = anom_disch)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#023E8A", "#E9A21B")) +  # colorblind-safe navy/amber
  labs(
    y = "Streamflow <br> (% difference from long-term mean)",
    x = "Year"
  ) +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(anom_disch$wyear) / 10),
      10 * ceiling(max(anom_disch$wyear) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(anom_disch$wyear) / 5),
      5 * ceiling(max(anom_disch$wyear) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)  # 5-yr minor ticks alongside 10-yr labeled major ticks
  ) +
  scale_y_symmetric() +
  geom_text(
    data = row_labels,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  facet_wrap(~local_site, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    axis.title.y = element_markdown(),  # renders the bold/italic markdown in the y title
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)  # extra right margin for the row label
  ) +
  # xlim pinned to the actual data range so the off-panel row label doesn't
  # stretch the x-axis; clip = "off" lets that label draw outside the panel
  coord_cartesian(xlim = x_range, clip = "off")

ggsave(g1,
       file = file.path(figures_dir, "streamflow_anom.png"),
       scale = 0.6, width = 10, height = 6
)


# ============================================================================
# FIGURE 2: Monthly streamflow anomaly at Green Lake 4 (May-Sep)
# ============================================================================
# Same anomaly logic as Figure 1, but computed per calendar month rather
# than per water year, and restricted to Green Lake 4.

anom_disch_month_gl4 <- infilled_df %>%
  filter(local_site == "gl4" & wyear < 2025) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(5:9)) %>%
  group_by(wyear, month) %>%
  summarise(total_disch = sum(discharge)) %>%
  ungroup() %>%
  group_by(month) %>%
  mutate(
    avgdisch = mean(total_disch, na.rm = TRUE),
    anom_disch = total_disch * 100 / avgdisch - 100,
    posneg = ifelse(anom_disch > 0, "pos", "neg") %>% factor(c("pos", "neg"))
  ) %>%
  ungroup() %>%
  mutate(month_name = month(month, label = TRUE)) %>%
  mutate(month_name = factor(month_name,
                             levels = c("May", "Jun", "Jul", "Aug", "Sep"), ordered = TRUE
  ))

x_range_month <- range(anom_disch_month_gl4$wyear, na.rm = TRUE)

# Row label on the rightmost panel of each row in the 3x2 facet grid
# (May/Jun/Jul on top, Aug/Sep on bottom -> Jul and Sep get the label)
row_labels_month <- data.frame(
  month_name = factor(c("Jul", "Sep"),
                      levels = c("May", "Jun", "Jul", "Aug", "Sep"), ordered = TRUE
  ),
  x = x_range_month[2] + diff(x_range_month) * 0.1,
  y = 0,
  label = "more~flow  %<->%  less~flow"
)

g2 <- ggplot(anom_disch_month_gl4, aes(x = wyear, y = anom_disch)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#023E8A", "#E9A21B")) +
  labs(
    y = "Streamflow <br> (% difference from long-term mean)",
    x = "Year"
  ) +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(anom_disch_month_gl4$wyear) / 10),
      10 * ceiling(max(anom_disch_month_gl4$wyear) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(anom_disch_month_gl4$wyear) / 5),
      5 * ceiling(max(anom_disch_month_gl4$wyear) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_symmetric() +
  geom_text(
    data = row_labels_month,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 3.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  facet_wrap(~month_name, nrow = 2, ncol = 3, axes = "all") +  # axes="all": every panel shows its own x-axis
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    axis.title.y = element_markdown(),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)
  ) +
  coord_cartesian(xlim = x_range_month, clip = "off")

ggsave(g2,
       file = file.path(figures_dir, "streamflow_anom_gl4_month.png"),
       scale = 0.7, width = 10, height = 7
)


# ============================================================================
# FIGURE 3: 2024 summer daily streamflow at Green Lake 4 vs. historical range
# ============================================================================
# Shows the current year's daily discharge as a black line, overlaid on a
# shaded 10th-90th percentile band built from all prior years' data for
# the same day-of-year.

summer_disch_last_yr <- infilled_df %>%
  filter(local_site == "gl4" & wyear < 2025) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(5:9)) %>%
  mutate(yday = lubridate::yday(date)) %>%
  filter(wyear == 2024) %>%
  left_join(
    infilled_df %>%
      filter(local_site == "gl4" & wyear < 2025) %>%
      mutate(yday = lubridate::yday(date), month = lubridate::month(date)) %>%
      filter(month %in% c(5:9)) %>%
      filter(wyear != 2024) %>%
      select(yday, discharge) %>%
      rename(oth = discharge)  # 'oth' = other (prior) years' discharge, for the ribbon
  )

hist_data_disch <- infilled_df %>%
  filter(local_site == "gl4" & wyear < 2025 & wyear != 2024)

fill_legend_disch <- paste0(
  "10th-90th percentile \n", "(",
  min(hist_data_disch$wyear), "-",
  max(hist_data_disch$wyear), ")"
)

g3 <- ggplot(summer_disch_last_yr %>% ungroup(), aes(x = date)) +
  geom_flquantiles(
    aes(y = oth),
    probs = c(0.10, 0.90), fill = "#4C72B0", alpha = 0.25, na.rm = TRUE
  ) +
  # Invisible single-row ribbon layer, present only to generate a legend
  # swatch for the historical band (geom_flquantiles doesn't support a
  # mapped fill aesthetic, so it can't create its own legend entry).
  geom_ribbon(
    data = summer_disch_last_yr %>% ungroup() %>% slice(1),
    aes(x = date, ymin = discharge, ymax = discharge, fill = fill_legend_disch),
    alpha = 0.25
  ) +
  geom_line(aes(y = discharge, color = as.character(wyear))) +
  scale_fill_manual(name = NULL, values = setNames("#4C72B0", fill_legend_disch)) +
  scale_color_manual(
    name = NULL,
    values = setNames("black", as.character(unique(summer_disch_last_yr$wyear)))
  ) +
  ylab("Streamflow (cubic meters)") +
  theme_hc() +
  labs(
    title = "Green Lake 4 Summer Streamflow",
    #subtitle = "Daily streamflow vs. historical range"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5)
  )

ggsave(g3,
       file = file.path(figures_dir, "streamflow_summer_gl4.png"),
       scale = 0.6, width = 8, height = 6
)

