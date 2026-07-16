# ============================================================================
# Niwot Ridge LTER Data Dashboard: Saddle Vegetation Growth Form Figures
# ============================================================================
#
# Purpose:
#   Reads Saddle grid point-quadrat vegetation survey data, classifies each
#   hit by growth form (forb, graminoid, shrub, moss, lichen), and produces:
#
#     1. sdl_veg_anom.jpg               Percent cover anomaly by growth form,
#                                        relative to each form's long-term mean
#     2. sdl_veg_abs.jpg                Absolute percent cover by growth form
#     3. sdl_veg_abs_horizontal.jpg     Same as (2), panels in a single row
#
# Input:
#   data/saddptqd.hh.data.csv     -- Saddle grid point-quadrat vegetation
#                                     survey, one row per point/hit
#                                     (EDI knb-lter-nwt.93)
#   data/pspecies_clean.js.data.csv -- species list with growth form
#                                     assignments (EDI knb-lter-nwt.91)
#
# Output:
#   Three JPG figures written to figures_dir (see file names above).
#
# TODO
#   - some USDA codes still lack a growth_habit assignment; check
#     unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)]) after any
#     data refresh and extend the manual assignment rules below if needed
#
# Original analysis: M. Oldfather, revised by Sarah C. Elmendorf,
#   10 Jan 2022, updated 14 Aug 2025
# Revised and updated for the Data Dashboard: Anne Marie Panetta, 2026
# ============================================================================

# -- SETUP --------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(svglite)
library(EDIutils)
library(lemon)
library(ggthemes)

options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# Anchor the working directory explicitly. 
# Update this if the project folder is ever moved or cloned elsewhere.
setwd("/Users/lter/Documents/GitHub/NWT LTER Data Dashboard/sdl_veg")


data_dir <- "data"
figures_dir <- "figures"

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
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  scope <- "knb-lter-nwt"  # Niwot scope
  
  # Data packages, with current citations (update as package versions change):
  #
  #   93 -- Saddle grid species composition
  #     Walker, M., H. Humphries, and Niwot Ridge LTER. 2025. Plant species
  #     composition data for Saddle grid, 1989 - ongoing. ver 10.
  #     Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/1427abaa317306bd0524046c0276b708
  #
  #   91 -- growth forms
  #     Smith, J., H. Humphries, M. Walker, and Niwot Ridge LTER. 2025.
  #     Plant species list for Niwot Ridge and Green Lakes Valley,
  #     1970 - ongoing. ver 3. Environmental Data Initiative.
  #     https://doi.org/10.6073/pasta/e96e1b8ce6d356154fb3a85eab96a106
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this loop.
  for (id in c("93", "91")) {
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

sdlcomp <- read.csv(file.path(data_dir, "saddptqd.hh.data.csv"), na.strings = na_vals)
sdlgf <- read.csv(file.path(data_dir, "pspecies_clean.js.data.csv"), na.strings = na_vals)


# -- CLEAN DATA -----------------------------------------------------------
# One growth-form error identified Jan 2026, not yet fixed on the portal:
# a lichen was mislabeled as litter for plot 1 in 2024.

sdlcomp <- mutate(
  sdlcomp,
  USDA_name = case_when(
    (USDA_name == "Litter" & plot == 1 & year == 2024) ~ "Lichen",
    TRUE ~ USDA_name
  ),
  USDA_code = case_when(
    (USDA_code == "2LTR" & plot == 1 & year == 2024) ~ "2LICHN",
    TRUE ~ USDA_code
  )
)

# Join in growth form, keeping only one row per USDA code (the species list
# can contain duplicate entries for a given code).
gf_join <- sdlgf %>%
  select(USDA_code, family, category, group, growth_habit) %>%
  distinct() %>%
  filter(!duplicated(USDA_code) & !duplicated(USDA_code), fromLast = TRUE)

sdlcomp <- sdlcomp %>%
  left_join(., gf_join)

# Spot-check which USDA names still lack a growth form after the join --
# these get manually assigned below.
unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit) | sdlcomp$growth_habit == ""])
unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])

# Manually assign growth form for entries not resolved by the species-list
# join (non-plant hits like rock/soil/scat, and a handful of plant genera
# not yet catalogued). Could be done more efficiently with case_when, but
# this ordered mutate() approach works fine for the current data size.

sdlcomp <- sdlcomp %>%
  mutate(
    growth_habit = ifelse((is.na(growth_habit) & grepl("Rock|Bare|Hole|soil", USDA_name)), "soil", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) & grepl("Scat|scat", USDA_name)), "scat", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) & grepl("Carex|Festuca", USDA_name)), "graminoid", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) &
      grepl(
        "^Polyg|^Arenaria|^Draba|^Erigeron|^Gentian|^Hieracium|^Noccaea|^Oreoxis|^Potentilla|^Ranunculus|^Saxifraga|^Silene|^Solidago|^Trifolium|^Stellaria|^Tetraneuris",
        USDA_name
      )), "forb", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) & grepl("Litter", USDA_name)), "litter", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) & grepl("Unknown", USDA_name)), "unknown", growth_habit),
    growth_habit = ifelse((is.na(growth_habit) & grepl("marker", USDA_name)), "marker", growth_habit)
  )

# Some USDA names still won't have a growth form assigned at this point --
# check the TODO note at the top of the file if this list is non-empty.
# unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])

# Plot 37 was sampled one year late in a couple of cases; treat those
# 1996 records as if collected in 1995 for consistency with other plots.

sdlcomp[sdlcomp$year == 1996, "year"] <- 1995

# Subset to one hit per point: the top hit if present, otherwise the
# bottom hit.
sdl_top <-
  sdlcomp %>%
  filter(hit_type == "bottom" | hit_type == "top") %>%
  arrange(year, plot, x, y, desc(hit_type)) %>%
  group_by(year, plot, x, y) %>%
  slice(1) %>%
  ungroup() %>%
  # add subshrubs to shrubs (Salix nivalis)
  mutate(growth_habit = ifelse(grepl("shrub", growth_habit), "shrub", growth_habit))

# -- SUMMARIZE COVER BY GROWTH FORM ----------------------------------------


# Count hits per growth form, per plot per year.
number_hits <- sdl_top %>%
  # add subshrubs to shrubs
  group_by(year, plot, growth_habit) %>%
  summarise(hits = dplyr::n(), .groups = "drop")

# Fill in explicit zeros for growth forms not hit in a given plot/year,
# so averages below aren't biased upward by only counting years/plots
# where a form was actually present.
data_with_zeros <- expand_grid(
  year = unique(number_hits$year),
  plot = unique(number_hits$plot),
  growth_habit = unique(number_hits$growth_habit)
)


number_hits <- left_join(data_with_zeros, number_hits)
number_hits[is.na(number_hits$hits), "hits"] <- 0


# Average across plots, by year and growth form.
number_hits_by_year <- number_hits %>%
  # add subshrubs to shrubs (Salix nivalis)
  group_by(year, growth_habit) %>%
  summarise(mean_cover = mean(hits))

# Quick look at all growth forms' trends before filtering to the ones used
# in the figures below.
ggplot(
  number_hits_by_year,
  aes(x = year, y = mean_cover, group = growth_habit, color = growth_habit)
) +
  geom_line() +
  facet_wrap(~growth_habit)

# -- DERIVE ANOMALIES -------------------------------------------------------

# Anomaly = each year's mean cover minus that growth form's long-term mean
# cover. Restricted to the five growth forms shown in the dashboard
# figures; "nonvascular"/"lichenous" are relabeled to the more familiar
# "moss"/"lichen".

anom_growth_form <- number_hits_by_year %>%
  group_by(growth_habit) %>%
  mutate(
    avgcov = mean(mean_cover, na.rm = TRUE),
    anom_cov = mean_cover - avgcov,
    posneg = ifelse(anom_cov > 0, "pos", "neg") %>%
      factor(c("pos", "neg")),
    growth_habit = ifelse(growth_habit == "nonvascular", "moss", growth_habit),
    growth_habit = ifelse(growth_habit == "lichenous", "lichen", growth_habit)
  ) %>%
  filter(growth_habit %in% c(
    "forb", "graminoid",
    "shrub", "moss", "lichen"
  )) %>%
  mutate(growth_habit = factor(growth_habit,
    levels = c(
      "forb", "graminoid",
      "shrub", "moss", "lichen"
    )
  ))

# ============================================================================
# FIGURE 1: Percent cover anomaly by growth form
# ============================================================================
# Each panel shows one growth form's cover anomaly (difference from its own
# long-term mean), colored by direction.

g1 <- ggplot(anom_growth_form, aes(x = year, y = anom_cov)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  labs(y = "Percent cover \n (Difference from long-term mean)", x = "Year", title = "Percent cover anomalies by growth form") +
  scale_x_continuous(
    breaks = seq(
      10 * floor(min(anom_growth_form$year) / 10),
      10 * ceiling(max(anom_growth_form$year) / 10),
      by = 10
    ),
    minor_breaks = seq(
      5 * floor(min(anom_growth_form$year) / 5),
      5 * ceiling(max(anom_growth_form$year) / 5),
      by = 5
    ),
    guide = guide_axis(minor.ticks = TRUE)  # 5-yr minor ticks alongside 10-yr labeled major ticks
  ) +
  theme_hc() +
  facet_wrap(~growth_habit, scales = "free_y") +
  theme(
    legend.position = "none",
    axis.minor.ticks.length = rel(0.5),
    plot.title = element_text(hjust = 0, margin = margin(b = 20)),
    plot.title.position = "plot"
  )

ggsave(g1,
       file = file.path(figures_dir, "sdl_veg_anom.jpg"),
       scale = 0.5, width = 15, height = 10, dpi = 600
)

# ============================================================================
# FIGURE 2: Absolute percent cover by growth form
# ============================================================================
# Companion figure to Figure 1, showing raw cover rather than anomalies, so
# growth forms can also be compared in absolute terms. Extra headroom is
# added above each panel so the highest bar isn't flush with the panel edge.

g2 <- anom_growth_form |>
  ggplot(aes(x = year, y = mean_cover, fill = growth_habit)) +
  geom_col(width = 0.75, alpha = 0.95, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.4) +
  scale_fill_manual(
    values = c(
      forb = "#D55EAA",
      graminoid = "#CFCF6B",
      shrub = "#006400",
      moss = "#8DAA3F",
      lichen = "#66C2A5"
    ),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23))) +  # ~5% white space above the tallest bar
  labs(y = "Percent cover", x = "Year", title= "Saddle vegetation cover by growth form") +
  theme_hc() +
  facet_wrap(~growth_habit, scales = "free_y") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = rel(0.95)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 8, b = 6, l = 0),
    plot.title = element_text(hjust = 0, margin = margin(t=10, b = 30)),
    plot.title.position = "plot"
  )

ggsave(g2,
       file = file.path(figures_dir, "sdl_veg_abs.jpg"),
       scale = 0.5, width = 15, height = 10, dpi = 600
)

# Alternate layout of Figure 2 for different display contexts (e.g. a
# single wide row for a dashboard banner, or a single tall column for a
# narrow sidebar).
ggsave(g2 + facet_wrap(~growth_habit, scales = "free_y", nrow = 1),
       file = file.path(figures_dir, "sdl_veg_abs_horizontal.jpg"),
       scale = 0.8, width = 10, height = 5, dpi = 600
)

