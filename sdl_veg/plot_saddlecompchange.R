# plot sdl spp comp by growth form
# SCE 10 Jan 2022
# updated 14 Jan 2025

# Credit for basic workflow for downloading/top subsetting/cleaning from
# M Oldfather script

# -- SETUP -----
# clean up enviro, read in needed libraries
rm(list = ls())
# library(readxl)
library(tidyverse)
library(svglite)
# library(vegan)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# only need to download once
download_data <- TRUE

# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data) {
  # download the data from EDI
  scope <- "knb-lter-nwt" # Niwot scope

  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c(
    "93", # saddle grid spp comp
    "91" # growth forms
  )) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")

    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")

    # download the data
    read_data_package_archive(packageID, path = "sdl_veg/data")
    print(read_data_package_citation(packageID))
    
    #[1] "Walker, M., H. Humphries, and Niwot Ridge LTER. 2025. Plant species composition data for Saddle grid, 1989 - ongoing. ver 10. Environmental Data Initiative. https://doi.org/10.6073/pasta/1427abaa317306bd0524046c0276b708. Accessed 2025-08-14."
    # "Smith, J., H. Humphries, M. Walker, and Niwot Ridge LTER. 2025. Plant species list for Niwot Ridge and Green Lakes Valley, 1970 - ongoing. ver 3. Environmental Data Initiative. https://doi.org/10.6073/pasta/e96e1b8ce6d356154fb3a85eab96a106. Accessed 2025-08-14."
  }

  # overwrites the manifests but don't really need them.
  for (fname in list.files("sdl_veg/data",
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = "sdl_veg/data/")
  }
}

# saddle comp by growth form ----------------------------------------------

sdlcomp <- read.csv("sdl_veg/data/saddptqd.hh.data.csv", na.strings = na_vals)
sdlgf <- read.csv("sdl_veg/data/pspecies_clean.js.data.csv", na.strings = na_vals)


# subset to usda only not duplicated
gf_join <- sdlgf %>%
  select(USDA_code, family, category, group, growth_habit) %>%
  distinct() %>%
  filter(!duplicated(USDA_code) & !duplicated(USDA_code), fromLast = TRUE)

sdlcomp <- sdlcomp %>%
  left_join(., gf_join)

unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit) | sdlcomp$growth_habit == ""])
unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])

# could more efficiently do with case_when but ordered mutate works ok for now
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

# some don't yet have growth form, assign
# unique (sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])

# change year = 1996 to 1995 for plot 37 (sampled one year late for a couple of plots)
sdlcomp[sdlcomp$year == 1996, "year"] <- 1995

# subset to only top hit (or bottom hit if no top hit present), by point
sdl_top <-
  sdlcomp %>%
  filter(hit_type == "bottom" | hit_type == "top") %>%
  arrange(year, plot, x, y, desc(hit_type)) %>%
  group_by(year, plot, x, y) %>%
  slice(1) %>%
  ungroup() %>%
  # add subshrubs to shrubs (Salix nivalis)
  mutate(growth_habit = ifelse(grepl("shrub", growth_habit), "shrub", growth_habit))


# aggregate data_file to calculate the number of hits
number_hits <- sdl_top %>%
  # add subshrubs to shrubs
  group_by(year, plot, growth_habit) %>%
  summarise(hits = dplyr::n())

# Fill in zeros for species/non-species not hit in a certain plot/year
data_with_zeros <- expand_grid(year = unique(number_hits$year), plot = unique(number_hits$plot), growth_habit = unique(number_hits$growth_habit))
data_with_zeros

# make abundance data frame
number_hits <- left_join(data_with_zeros, number_hits)
number_hits[is.na(number_hits$hits), "hits"] <- 0


# aggregate w the zeros
number_hits_by_year <- number_hits %>%
  # add subshrubs to shrubs (Salix nivalis)
  group_by(year, growth_habit) %>%
  summarise(mean_cover = mean(hits))

ggplot(
  number_hits_by_year,
  aes(x = year, y = mean_cover, group = growth_habit, color = growth_habit)
) +
  geom_line() +
  facet_wrap(~growth_habit)


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

# plot
g1 <- ggplot(anom_growth_form, aes(x = year, y = anom_cov)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  labs(y = "Percent cover \n difference from long-term mean", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(more %<->% less))) +
  theme_hc() +
  facet_wrap(~growth_habit, scales = "free_y") +
  theme(legend.position = "none")

ggsave(g1,
  file = file.path("sdl_veg", "figures", "sdl_veg.jpg"),
  scale = 0.8, width = 10, height = 6, dpi = 600
)

