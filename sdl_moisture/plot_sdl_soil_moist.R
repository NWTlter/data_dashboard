# Script to download and plot Saddle (nwt package id 405) soil moisture
# Code developed by Sarah Elmendorf. 
# Updated for Data Dashbaord by Anne Marie Panetta 2026.
# temporal patterns

library(tidyverse)
library(ggthemes)
library(lemon)
library(officer)
library(magrittr)
library(viridisLite) #for color-blind accessibility
library(EDIutils)

rm(list = ls())

# Define directory paths for cross-platform compatibility
# Go up one level to data_dashboard, then into sdl_moisture
parent_dir <- dirname(getwd())  # Go from c1_d1_sdl_temp_ppt to data_dashboard
data_dir <- file.path(parent_dir, "sdl_moisture", "data")
figures_dir <- file.path(parent_dir, "sdl_moisture", "figures")

# Create directories if they don't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Verify paths
cat("\n=== Directory Setup ===\n")
cat("Current working dir:", getwd(), "\n")
cat("Parent dir:", parent_dir, "\n")
cat("Data dir:", data_dir, "\n")
cat("  Exists:", dir.exists(data_dir), "\n")
cat("Figures dir:", figures_dir, "\n")
cat("  Exists:", dir.exists(figures_dir), "\n")
cat("========================\n\n")

# only need to download once
download_data <- FALSE

# User-specified notes for PowerPoint slides
# slide g1 - no dynamic text edit directly
slide_notes_g1 <- "This figure below illustrates how summer soil moisture (June-August) deviates from average across the Niwot Ridge alpine tundra soil moisture record. Data are from an alpine dry meadow located on the Saddle of Niwot Ridge."

# Define text for slide notes g2
# dynamic text, the years will be inserted later
slide_notes_g2_text <- list(
  "Current year soil moisture compared to historical patterns. The ribbon shows the 10th-90th percentile range of historical data (",
  # min(df$year), "-", max(df$year)-1,
  "), providing context for evaluating whether conditions in ",
  # max(df$year),
  " are within normal ranges."
)

citation_text <- "Morse, J. and M. Losleben. 2025. Climate data for saddle data loggers (CR23X and CR1000), 2000 - ongoing, daily. ver 10. Environmental Data Initiative. https://doi.org/10.6073/pasta/b01aea637f7608f0a1b2895ae474d571. Accessed 2025-08-11."

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
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }

  # update the below so you remember to cite it correctly
  # "Morse, J. and M. Losleben. 2025. Climate data for saddle data loggers (CR23X and CR1000), 2000 - ongoing, daily. ver 10. Environmental Data Initiative. https://doi.org/10.6073/pasta/b01aea637f7608f0a1b2895ae474d571. Accessed 2025-08-11."
  # overwrites the manifests but don't really need them.
  for (fname in list.files(data_dir,
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}

# read and munge ppt and temp --------------------------------------------------

df <- read.csv(file.path(data_dir, "sdlcr23x-cr1000.daily.ml.data.csv"),
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
  mutate(tmin = (soiltemp_5cm_min + 2) / (5 + 2)) %>%
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
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(y = "Soil moisture anomaly (%) \n Jun-Aug", x = "Year") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(wetter %<->% drier))) +
  theme_hc() +
  theme(legend.position = "none")
g1

ggsave(g1,
  file = file.path(figures_dir, "sdl_soil_moist_anom.png"),
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
    .groups = "drop"
  )

fill_legend <- paste0(
  "10th-90th percentile \n", "(",
  min(df$year), "-",
  max(df$year) - 1, ")"
)

g2 <- ggplot(summer_moist_last_yr %>%
  ungroup(), aes(x = yday)) +
  geom_line(aes(
    y = soilmoist_5cm_avg,
    color = as.character(year)
  )) +
  geom_ribbon(
    data = byday,
    aes(
      x = yday, ymin = q10, ymax = q90,
      fill = fill_legend
    ),
    alpha = 0.3
  ) +
  scale_x_continuous(
    breaks = c(121, 152, 182, 213, 244), # May 1, June 1, July 1, Aug 1, Sep 1 (approx)
    labels = c("May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1")
  ) +
  scale_color_manual(
    name = NULL,
    values = setNames("black", as.character(unique(summer_moist_last_yr$year)))
  ) +
  scale_fill_manual(
    name = NULL,
    values = setNames("pink", fill_legend)
  ) +
  ylab("Soil moisture (%)") +
  xlab(NULL) +
  theme_hc()
g2

ggsave(g2,
  file = file.path(figures_dir, "sdl_soil_moist_last_year.png"),
  scale = 0.5, width = 8, height = 6
)


# Create first PowerPoint presentation for g1
ppt1 <- read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  # ph_with(value = "SDL Soil Moisture Analysis", location = ph_location_type(type = "title")) %>%
  ph_with(
    value = external_img(file.path(figures_dir, "sdl_soil_moist_anom.png")),
    location = ph_location(left = 1, top = 0.5, width = 8, height = 5.5)
  ) %>%
  ph_with(
    value = fpar(
      ftext("Saddle (SDL) Data: ", prop = fp_text(bold = TRUE, font.size = 12, font.family = "sans")),
      ftext(citation_text, prop = fp_text(font.size = 12, font.family = "sans"))
    ),
    location = ph_location(left = 0.5, top = 6.5, width = 9, height = 1)
  ) %>%
  set_notes(value = slide_notes_g1, location = notes_location_type(type = "body"))

# Define slide_notes_g2 after df is loaded
slide_notes_g2 <- paste0(slide_notes_g2_text[[1]], min(df$year), "-", max(df$year) - 1, slide_notes_g2_text[[2]], max(df$year), slide_notes_g2_text[[3]])

# Create second PowerPoint presentation for g2
ppt2 <- read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  # ph_with(value = "Current Year Soil Moisture Patterns", location = ph_location_type(type = "title")) %>%
  ph_with(
    value = external_img(file.path(figures_dir, "sdl_soil_moist_last_year.png")),
    location = ph_location(left = 1, top = 0.5, width = 8, height = 5.5)
  ) %>%
  ph_with(
    value = fpar(
      ftext("Saddle (SDL) Data: ", prop = fp_text(bold = TRUE, font.size = 12, font.family = "sans")),
      ftext(citation_text, prop = fp_text(font.size = 12, font.family = "sans"))
    ),
    location = ph_location(left = 0.5, top = 6.5, width = 9, height = 1)
  ) %>%
  set_notes(value = slide_notes_g2, location = notes_location_type(type = "body"))

# Save both PowerPoint presentations
print(ppt1, target = file.path(figures_dir, "sdl_soil_moist_anom.pptx"))
print(ppt2, target = file.path(figures_dir, "sdl_soil_moist_last_year.pptx"))

cat("PowerPoint presentations saved to:\n")
cat("  -", file.path(figures_dir, "sdl_soil_moist_anom.pptx"), "\n")
cat("  -", file.path(figures_dir, "sdl_soil_moist_last_year.pptx"), "\n")
