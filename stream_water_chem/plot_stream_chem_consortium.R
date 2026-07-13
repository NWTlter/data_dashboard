# Code to plot concentrations
# Code developed by Sarah Elmendorf (2023)
# Code updated by Anne Marie Panetta (2026)
# see nwt-8 long-term trends for flow-normalized workflow if desired

# setup#########
# clean up enviro, read in needed libraries
rm(list = ls())
# library(readxl)
library(tidyverse)
# library(vegan)
options(stringsAsFactors = F)
theme_set(theme_bw())

####New from Anne Marie to download data from EDI####
# only need to download once
download_data <- TRUE

data_dir <- "data"
figures_dir <- "figures"

# Create figures directory if it doesn't exist
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# download data -----------------------------------------------------------
# Note: if you have already downloaded SOME data, the read_data_package_archive
# function will not work as it does not want to overwrite. Clear your /data
# directories and then return to run the following code.
if (download_data) {
  # download the data from EDI
  # discharge inputs used are:
  # knb-lter-nwt.102 (albion)
  # knb-lter-nwt.105 (gl4)
  # knb-lter-nwt.111 (martinelli)
  # knb-lter-nwt.74 (saddle)
  
  # no discharge but if you just want trends
  # there are some others in just select yrs, etc.
  # 163 is gr5rock glacier; 162 is watershed flume; 213 is soddie, et
  # 109 is gl5#rg
  
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  scope <- "knb-lter-nwt" # Niwot scope
  
  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this
  for (id in c(
    "103", "108", "160",
    "112", "163", "109"
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
# Albion: 
# GL4
# Mart:
# GL5:
# GrrG:
  
for (fname in list.files(data_dir,
                         pattern = "knb-lter.*zip",
                         full.names = TRUE
)) {
  unzip(zipfile = fname, exdir = data_dir)
}
}

# read in data --------------------------------------------------

# all missing values that have to are MCAR
na_vals <- c("", "NA", "NaN", "NP", "DNS", "NSS", "EQCL", "QNS", "NV", "dns")

data_file <- list()
data_file[["albisolu"]] <- read.csv(file.path(data_dir, "albisolu.nc.data.csv"),
                                    na = na_vals)
data_file[["gre4solu"]] <- read.csv(file.path(data_dir, "gre4solu.nc.data.csv"),
                                    na = na_vals)
data_file[["saddsolu"]] <- read.csv(file.path(data_dir, "saddsolu.nc.data.csv"),
                                    na = na_vals)
data_file[["martsolu"]] <- read.csv(file.path(data_dir, "martsolu.nc.data.csv"),
                                    na = na_vals)
data_file[["grrgsolu"]] <- read.csv(file.path(data_dir, "grrgsolu.nc.data.csv"),
                                    na = na_vals)
data_file[["gre5solu"]] <- read.csv(file.path(data_dir, "gre5solu.nc.data.csv"),
                                    na = na_vals)


datafile <- plyr::rbind.fill(data_file) %>% data.frame(.)

num_vars <- names(datafile)[6:44]

# define replace function for non-numeric values
repl.f <- function(x) ifelse(grepl("<|trace|u", x), 0, x)

# replace all trace, u, and <mdl values with 0 - this may
# not be the best thing to do but seems better than leaving them out?
datafile <- datafile %>% dplyr::mutate_at(vars(num_vars), .funs = repl.f)


# do by may, june (melt phase)
# july, aug sept

# nitrate, sulfateNO3./SO4..

# remove high values >3sd greater than that months lt-average
# as there seem to be some big outliers in there
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

# make sulfate anoms by month for late season sulfur
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


# make all into anomalies, by month
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

g1 <- ggplot(anom_chem_month, aes(x = year, y = anom_S)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#99000d", "#034e7b")) +
  labs(y = "September sulfate anomaly (%)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(
    trans = I,
    breaks = NULL, name = expression(more %<->% less)
  )) +
  theme_hc() +
  facet_wrap(~local_site) +
  theme(legend.position = "none")


ggsave(g1,
  file = "plots/sulfate_anom_sept.png",
  scale = 0.5, width = 8, height = 6
)


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

# make all into anomalies, by month or not
anom_chem_summer_N <- chem_by_summer %>%
  filter(local_site %in% c("ALBION", "GREEN LAKE 4", "GREEN LAKE 5", "MARTINELLI")) %>%
  # group_by(local_site, year) %>%
  # summarise(mean_N = mean (nitrate, na.rm = TRUE, .groups = 'drop'))%>%
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

g2 <- ggplot(anom_chem_summer_N, aes(x = year, y = anom_N)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#99000d", "#034e7b")) +
  labs(y = "Summer Nitrate anomaly (%)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(
    trans = I,
    breaks = NULL, name = expression(more %<->% less)
  )) +
  theme_hc() +
  facet_wrap(~local_site) +
  theme(legend.position = "none")


ggsave(g2,
  file = "plots/nitrate_anom_summer.png",
  scale = 0.5, width = 8, height = 6
)
