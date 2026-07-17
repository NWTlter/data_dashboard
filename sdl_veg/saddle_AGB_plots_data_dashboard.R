# ============================================================================
# Niwot Ridge LTER Data Dashboard: Saddle ANPP (AGB) Anomaly Figure
# ============================================================================
#
# Purpose:
#   Reads Saddle grid aboveground net primary productivity (ANPP) data and
#   produces:
#
#     1. sdl_aboveground_biomass_anom.png   Aboveground biomass (ANPP) anomaly, 
#                                           relative to the long-term mean
#
# Input:
#   data/saddanpp.hh.data.csv  -- Saddle grid ANPP (EDI knb-lter-nwt.16)
#
# Output:
#   One PNG figure written to figures_dir (see file name above).
#
#
# Credits: 
#   Sarah C. Elmendorf (2022) developed this script based on a download/
#   subsetting workflow by M. Oldfather (nwt_8-renewal/get_saddle_sp_CNM.R)
#   Anne Marie Panetta (2026) updated and revised this script for the 
#   data dashboard 
# ============================================================================

# -- SETUP --------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(EDIutils)  # tools for interacting with EDI's data package API
library(ggthemes)  # theme_hc()
library(ggpmisc)   # symmetric_limits() for scale_y_continuous

options(stringsAsFactors = FALSE)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# Anchor working directory to this script's location. Note: this only
# works in an interactive RStudio session with the script saved and the
# active/focused tab -- if you hit "cannot open file" errors downstream,
# check getwd() and consider hardcoding this path instead.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
  
  # 16 -- Saddle grid ANPP, with current citation (update as package
  # version changes):
  #   Walker, M., J. Smith, H. Humphries, and Niwot Ridge LTER. 2025.
  #   Aboveground net primary productivity data for Saddle grid,
  #   1992 - ongoing. ver 10. Environmental Data Initiative.
  #   https://doi.org/10.6073/pasta/5d013d5908bf787b7aabd5ba08ac71f5
  for (id in c("16")) {
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    packageID <- paste(scope, id, revision, sep = ".")
    
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID)) # confirm you're citing the version actually used
  }
  
  for (fname in list.files(data_dir, pattern = "knb-lter.*zip", full.names = TRUE)) {
    unzip(zipfile = fname, exdir = data_dir)
  }
  
  # After the first download, list what actually landed in data_dir and
  # update the file name in READ DATA below to match, if needed.
  print(list.files(data_dir, pattern = "\\.csv$", ignore.case = TRUE))
}


# -- READ DATA ------------------------------------------------------------
# TODO: confirm the actual file name in data_dir once downloaded (see
# list.files() above) -- this is a best-guess name based on the package
# title and other dashboard scripts' naming conventions.

sdlprod <- read.csv(file.path(data_dir, "saddgrid_npp.hh.data.csv"), na.strings = na_vals)


# ============================================================================
# FIGURE: Aboveground biomass (ANPP) anomaly
# ============================================================================
# Year-level ANPP anomaly, expressed as percent difference from the
# long-term mean. Grid points are averaged first (across the n=2 subsamples
# collected per point), then years are averaged across grid points.

npp_by_year <- sdlprod %>%
  group_by(year, grid_pt) %>%
  summarise(NPP = mean(NPP)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(NPP = mean(NPP)) %>%
  ungroup() %>%
  mutate(
    avgNPP = mean(NPP, na.rm = TRUE),
    anom_NPP = (NPP * 100 / avgNPP) - 100,
    posneg = ifelse(anom_NPP > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

x_range_npp <- range(npp_by_year$year, na.rm = TRUE)

row_label_npp <- data.frame(
  x = x_range_npp[2] + diff(x_range_npp) * 0.1,
  y = 0,
  label = "more~biomass  %<->%  less~biomass"
)

g1 <- ggplot(npp_by_year, aes(x = year, y = anom_NPP)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  labs(y = "Aboveground biomass \n (% difference from long-term mean)", x = "Year") +
  scale_x_continuous(
    breaks = seq(
      5 * floor(min(npp_by_year$year) / 5),
      5 * ceiling(max(npp_by_year$year) / 5),
      by = 5
    )
  ) +
  scale_y_continuous(limits = symmetric_limits) +
  geom_text(
    data = row_label_npp,
    aes(x = x, y = y, label = label),
    angle = -90, hjust = 0.5, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE, parse = TRUE
  ) +
  theme_hc() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)  # extra right margin for the row label
  ) +
  coord_cartesian(xlim = x_range_npp, clip = "off")

ggsave(g1,
       file = file.path(figures_dir, "sdl_aboveground_biomass_anom.png"),
       scale = 0.8, width = 10, height = 6
)