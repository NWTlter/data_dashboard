# Use interval censored model to more precisely estimate snowmelt dates
# from weekly/biweekly/sometimes more snow surveys
# SCE 11 July 2025
# note you will need 4 cores and some time to run this code

# -- SETUP ------
library(tidyverse) # munging
#library(lme4) # space_time of NPP vs snowmelt model
#library(emmeans) # not yet used but should keep to test mixed models eventually
#library(lmerTest) # not yet used but should keep to test mixed models eventually
library(brms) # gap fill snowmelt dates
library(EDIutils) # Handy tools for interacting with EDI's API

# only need to download once
download_data <- TRUE

#slow to run, might want to run only once
run_bayesian_models <-FALSE

# functions -----------------------------------------------------------

# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data) {
  # download the data from EDI
  # 31 is sdl snow
  
  data_dir <- file.path("sdl_snow", "data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  scope <- "knb-lter-nwt" # Niwot scope
  
  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c("31")) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    
    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")
    
    # download the data
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }
  
  # update the below so you remember to cite it correctly
  # [1] "Walker, D., J. Morse, and Niwot Ridge LTER. 2024. Snow depth data for Saddle grid, 1992 - ongoing. ver 21. Environmental Data Initiative. https://doi.org/10.6073/pasta/54fd300090a4859ad3083805e98bc823. Accessed 2025-07-11."
  # overwrites the manifests but don't really need them.
  for (fname in list.files(data_dir,
                           pattern = "knb-lter.*zip",
                           full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = "sdl_snow/data/")
  }
}


#Read in data####-----------------------------------------------------------------

# Bring in snow data
snow <- read_csv(file.path( data_dir, "saddsnow.dw.data.csv"))

# Seems as though they skipped one meas in a few years, but before/after was 0
#in that grid so assume it's still 0

snow <- snow %>%
  mutate(
    mean_depth = ifelse(
      date == "2000-06-09" & point_ID == 301,
      0, mean_depth
    ),
    mean_depth = ifelse(
      date == "1998-01-09" & point_ID == 10,
      0, mean_depth
    ),
    mean_depth = ifelse(
      date == "1998-05-18" & point_ID == 401,
      0, mean_depth
    ),
    date = ymd(date),
    year = year (date),
    month = month(date),
    doy = yday (date),
    season = ifelse(doy > 240, "fall", "spring"),
    # At the moment we're not using depth as an absolute
    # so we can ignore the "deeper than pole" uncertainties
    mean_depth = as.numeric(gsub('\\+', '', mean_depth))
  )



#Calculate melt dates####------------------------------------------------------------
last_snow_covered <- snow %>%
  filter(mean_depth > 0 & !is.na(mean_depth) & season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(last_snow = max(doy), .groups = "drop")

first_melt <- snow %>%
  full_join(., last_snow_covered) %>%
  filter(doy > last_snow & season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(snowmelt = min(doy), .groups = "drop")

# There are some strong outliers - It is unclear
# whether these points really had no snow then or
# the sampling skipped them. 
# It is also not clear what the biological meaning 
# of very early melt dates are; 
# these points did receive snow after the melt date,
# but this snow did not stick around long enough to be sampled
hist (first_melt$snowmelt)

# Some are 0 by the time sampling people arrive in a year,
# add those as right censored at 0.
# It is probably not actually the case that there
# was never snow. It is more likely that snow was intermittent, and samplers
# never arrived at a time when snow had not blown off.

no_snow <- snow %>%
  filter(season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(snowtot = sum(mean_depth), .groups = "drop") %>%
  filter(snowtot == 0) %>%
  left_join(snow %>%
              group_by(year, point_ID) %>%
              summarize(start_samp = min(doy)))

# Last sample not snowfree
no_snowfree <- snow %>%
  filter(season == "spring") %>%
  group_by(year, point_ID) %>%
  summarize(doy = max(doy)) %>%
  left_join(., snow) %>%
  filter(mean_depth != 0) %>%
  rename(last_snow_samp = doy)


# join
snow_cens <- last_snow_covered %>%
  full_join(., first_melt) %>%
  full_join(., (snow %>%
                  select(point_ID, year) %>%
                  distinct())) %>%
  full_join(., no_snow) %>%
  full_join(., no_snowfree) %>%
  rename(grid_pt = point_ID) %>%
  mutate(
    snowmelt =
      dplyr::if_else(!is.na(start_samp), start_samp, snowmelt)
  ) %>%
  mutate(
    last_snow =
      dplyr::if_else(!is.na(last_snow_samp), last_snow_samp, last_snow)
  ) %>%
  filter(year > 1992) 


# Another option is to reset the very early days to sometime
# between 75 and 105 (this is something like mid-march - mid april).
# This is the approach we will take here:
snow_cens_2 <- last_snow_covered %>%
  full_join(., first_melt) %>%
  full_join(., (snow %>%
                  select(point_ID, year) %>%
                  distinct())) %>%
  full_join(., no_snow) %>%
  full_join(., no_snowfree) %>%
  rename(grid_pt = point_ID) %>%
  mutate(
    snowmelt =
      dplyr::if_else(!is.na(start_samp), start_samp, snowmelt)
  ) %>%
  mutate(
    last_snow =
      dplyr::if_else(!is.na(last_snow_samp), last_snow_samp, last_snow)
  ) %>%
  mutate (snowmelt_adj = ifelse(snowmelt<105&!is.na(snowmelt), 105, snowmelt),
          last_snow_adj = if_else((last_snow<75 &!is.na(last_snow)),
                                  75, last_snow)
  )%>%
  filter(year > 1992)


if (run_bayesian_models){

# The censoring variable
# (named "censored" in this example) should contain the values 'left', 'none',
#' right', and 'interval' (or equivalently -1, 0, 1, and 2).
#'
#' Because I can never remember which is which - 
#' #Left censoring is when the event of interest has already occurred before
# enrollment. This is very rarely encountered. Truncation is deliberate and
# due to study design.

# center to make modeling happy
mn <- mean(c(snow_cens$snowmelt, snow_cens$last_snow), na.rm = TRUE)
sd <- sd(c(snow_cens$snowmelt, snow_cens$last_snow), na.rm = TRUE)

snow_cens <- snow_cens %>%
  mutate(
    censored =
      case_when(
        !is.na(last_snow) & !is.na(snowmelt) ~ "interval", # has both dates
        is.na(last_snow) & !is.na(snowmelt) ~ "left", # first visit has 0 snow
        is.na(snowmelt) & !is.na(last_snow) ~ "right", # last visit still has snow
        TRUE ~ "NA"
      ),
    year = factor(year),
    grid_pt = factor(grid_pt)
  ) %>%
  mutate(
    last_snow_c = (last_snow - mn) / sd,
    snowmelt_c = (snowmelt - mn) / sd
  ) %>%
  # for how brms wants the censored data coded
  mutate(y = ifelse(censored == "left", snowmelt_c, last_snow_c)) %>%
  # For technical reasons, y2 needs to be specified for all observations
  # even if they are unused, for left and none, for instance.
  # In the future brms will be better at selectively checking whether variables
  # are actually used in a specific case but this is hard. So for now:
  # Any variable that you use in your data should be non NA for all observation
  # you want to include, even if some of the values are actually unused in the model.
  # https://github.com/paul-buerkner/brms/issues/799
  mutate(snowmelt_c = ifelse(censored == "right", 999, snowmelt_c))

# run model, takes a while, should set up to use more cores at some point
fit_cens <- brm(y | cens(censored, snowmelt_c) ~ 1 + (1 | grid_pt) + (1 | year),
                data = snow_cens, iter = 20000, cores = 4
)

saveRDS (fit_cens, 'sdl_snow/big_bayesian_files/fit_cens_snow_cens.rds' ) 
#I think we need to update this using file.path, but I don't want to mess this up.


# Center to make modeling happy
mn_2 <- mean(c(snow_cens_2$snowmelt_adj, snow_cens_2$last_snow_adj), na.rm = TRUE)
sd_2 <- sd(c(snow_cens_2$snowmelt_adj, snow_cens_2$last_snow_adj), na.rm = TRUE)

snow_cens_2 <- snow_cens_2 %>%
  mutate(
    censored =
      case_when(
        !is.na(last_snow_adj) & !is.na(snowmelt_adj) ~ "interval", # has both dates
        is.na(last_snow_adj) & !is.na(snowmelt_adj) ~ "left", # first visit has 0 snow
        is.na(snowmelt_adj) & !is.na(last_snow_adj) ~ "right", # last visit still has snow
        TRUE ~ "NA"
      ),
    year = factor(year),
    grid_pt = factor(grid_pt)
  ) %>%
  mutate(
    last_snow_c = (last_snow_adj - mn_2) / sd_2,
    snowmelt_c = (snowmelt_adj - mn_2) / sd_2
  ) %>%
  # for how brms wants the censored data coded
  mutate(y = ifelse(censored == "left", snowmelt_c, last_snow_c)) %>%
  # For technical reasons, y2 needs to be specified for all observations
  # even if they are unused, for left and none, for instance.
  # In the future, brms will be better at selectively checking whether variables
  # are actually used in a specific case but this is hard. So for now:
  # Any variable that you use in your data should be non NA for all observations
  # you want to include, even if some of the values are actually unused in the model.
  # https://github.com/paul-buerkner/brms/issues/799
  mutate(snowmelt_c = ifelse(censored == "right", 999, snowmelt_c))

# Run model (this takes a while; we should set this up to use more cores at some point)
fit_cens_2 <- brm(y | cens(censored, snowmelt_c) ~ 1 + (1 | grid_pt) + (1 | year),
                  data = snow_cens_2, iter = 20000, cores = 4)

saveRDS (fit_cens_2, 'sdl_snow/big_bayesian_files/fit_cens_snow_cens_2.rds' )

}

# Read in the bayesian model results and continue ------------------------------
fit_cens_2 <-readRDS('sdl_snow/big_bayesian_files/fit_cens_snow_cens_2.rds')
#The above also needs to be updated to use file.path(), I think

# Predictions
pp_2 <- predict(fit_cens_2, newdata = snow_cens_2)

pp_2 <- pp_2 %>%
  data.frame() %>%
  mutate(
    Estimate = (Estimate * sd_2) + mn_2,
    Q2.5 = (Q2.5 * sd_2) + mn_2,
    Q97.5 = (Q97.5 * sd_2) + mn_2
  )

pp_2 <- cbind(pp_2, snow_cens_2)

ggplot(pp_2, aes(x = last_snow_adj, y = Estimate)) +
  geom_point() +
  geom_abline()

# This seems to work ok, but some predicted values fall outside the real range.
# Basically, it is hard to predict the very windblown sites.
ggplot(pp_2, aes(
  x = Estimate, ymin = last_snow_adj, ymax = snowmelt_adj,
  color = year
)) +
  geom_linerange() +
  geom_abline()

# Solution - 
# Use estimate at face value if in bounds, otherwise, take bounds

pp_2 <- pp_2 %>%
  mutate(
    snowmelt_est =
      case_when(
        !is.na(snowmelt_adj) & !is.na(last_snow_adj) & Estimate > last_snow_adj & Estimate < snowmelt_adj ~ Estimate,
        !is.na(snowmelt_adj) & is.na(last_snow_adj) & Estimate < snowmelt_adj ~ Estimate,
        is.na(snowmelt_adj) & !is.na(last_snow_adj) & Estimate > last_snow_adj ~ Estimate,
        !is.na(snowmelt_adj) & !is.na(last_snow_adj) & (Estimate > last_snow_adj | Estimate < snowmelt_adj) ~ (snowmelt_adj + last_snow_adj) / 2,
        !is.na(snowmelt_adj) & is.na(last_snow_adj) & Estimate > snowmelt_adj ~ snowmelt_adj,
        is.na(snowmelt_adj) & !is.na(last_snow_adj) & Estimate < last_snow_adj ~ last_snow_adj
      )
  )

write.csv(pp_2 %>%
            select(year, grid_pt,last_snow_adj, snowmelt_adj, snowmelt_est, Estimate, Est.Error, Q2.5, Q97.5),
          file = 'sdl_snow/data_deriv/snowmelt_est_adj.csv',
          row.names = FALSE)
# We may need to update to file path here too. 