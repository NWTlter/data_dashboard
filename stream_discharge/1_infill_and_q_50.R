# code to infill and play around with discharge phenology from NWT
# sce 20 Jan 2021.

# todo####
# consider keeping the records w/
# ""estimate based on partial records only"

# consider weighing analyses based on amt
# of infilled data or removing largely infilled years

# updated code to read from latest version on EDI.


# setup###################################################################
# libraries
library(ggplot2)
library(tidyverse)
library(mtsdi) # for imputing

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# try a few different sp_controls, they don't in the end matter that much
# here I use 7, but tried 3,5,9 just for a viz check
sp_control <- 7

# will makes plots dump out to a file if false, which is kind of helful
plot_in_code <- FALSE


# read data###########################################################
# read in all data - latest versions from nwt server
# at some point update this to read off EDI
all_files <- c(
  "albdisch.nc.data.csv",
  "gl4disch.nc.data.csv",
  "mardisch.nc.data.csv",
  "saddisch.nc.data.csv"
)

df <- list()
for (file in all_files) {
  df[[file]] <- readr::read_csv(paste0("http://niwot.colorado.edu/data_csvs/", file),
    na = "NaN", guess_max = 1000000
  )
}

df <- df %>%
  data.table::rbindlist(., use.names = TRUE, idcol = "file", fill = TRUE)

### infill martinelli and saddle###############################################

# first add in notes from mart in years where there are real notes
# on when it went dry
# check against mart
# sce some of these are in there but some gaps could fill
df <- df %>%
  # 1998 Record ends 09/24/98 when channel dry and logger removed.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "1998-09-24" & date < "1999-01-01",
  0, discharge
  )) %>%
  # 1999 Record ends when channel dry and logger removed (13:00 MST 14/09/99).
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "1999-09-14" & date < "2000-01-01",
  0, discharge
  )) %>%
  # 2012: The channel was dry after 9/08/12.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2012-08-09" & date < "2013-01-01",
  0, discharge
  )) %>%
  # 2015: The channel was dry after 10/09/15.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2015-10-15" & date < "2016-01-01",
  0, discharge
  )) %>%
  # 2015: The channel was dry after 10/09/15.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2015-10-15" & date < "2016-01-01",
  0, discharge
  )) %>%
  # 2017: The channel was dry after 21 September 2017 17:00 MST.
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2017-09-17" & date < "2018-01-01",
  0, discharge
  )) %>%
  # 2018: Channel was dry after 16 August 2018 13:25 MST
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2018-08-16" & date < "2019-01-01",
  0, discharge
  )) %>%
  # 2019 Channel was dry after 21 August 2019
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
    date >= "2019-08-21" & date < "2020-01-01",
  0, discharge
  )) %>%
# 2020 Channel was dry after 24 Sept 2020
mutate(discharge = ifelse(is.na(discharge) & local_site == "mar" &
  date >= "2020-09-24" & date < "2021-01-01",
0, discharge
))



# Mart and saddle are seasonal, so first step is to infill
# the doys where there is NEVER flow w 0s
# calculate flow season
flow_season <- df %>%
  mutate(yday = lubridate::yday(date)) %>%
  mutate(discharge = ifelse(!is.na(notes), NA, discharge)) %>%
  filter(!is.na(discharge) & discharge > 0) %>%
  group_by(local_site) %>%
  summarise(start = min(yday), end = max(yday))

# fill mart and saddle consistent off season with 0
infilled_df <- df %>%
  mutate(
    yday = lubridate::yday(date),
    year = lubridate::year(date)
  ) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
  left_join(flow_season) %>%
  # do not use any previously infilled data
  mutate(discharge = ifelse(!is.na(notes) & !notes %in% c(
    "estimate based on partial records only",
    "flow data estimated from observations",
    "flow data estimated from intermittent observation",
    "flow data estimated from intermittent observations",
    "flow data estimated from field observations"
  ), NA, discharge)) %>%
  mutate(discharge = ifelse(yday < start, 0, discharge)) %>%
  mutate(discharge = ifelse(yday > end, 0, discharge)) # %>%
# left_join(., gap_by_wy)

# step 2 of mar/saddle is to use mnimmput with
# the other site as a predictor, as well as all other
# years at the same site as a predictor
# but truncating out the winter season

# saddle first

df_wide <- infilled_df %>%
  filter(yday > min(flow_season$start[flow_season$local_site %in% c("sdl", "mar")]) &
    yday < max(flow_season$end[flow_season$local_site %in% c("sdl", "mar")])) %>%
  filter(local_site %in% c("mar", "sdl")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
  arrange(date)

# add in the sdl for all years
for (yr in unique(df_wide$year)) {
  newdat <- df %>%
    filter(local_site == "sdl" & lubridate::year(date) == yr) %>%
    mutate(yday = lubridate::yday(date))
  if (yr == 2013) {
    newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove flood as seasonal predictor
  }
  # skip as predicators any years that are all na
  if (!all(is.na(newdat$discharge))) {
    newdat[[paste0("sdl_", yr)]] <- newdat[["discharge"]]
    df_wide <- df_wide %>%
      left_join(., newdat %>%
        select(
          -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
          -notes
        ))
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$year)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$year == yr] <- NA
  }
}

df_wide <- df_wide %>%
  # remove flood
  mutate(sdl = ifelse(year == 2013 & yday %in% c(253:259), NA, sdl)) %>%
  arrange(date)

# make predictions
# multivariate EM imputation w ts, spline w 7df

pred_sdl <- mnimput(df_wide %>%
  select(-date, -year, -yday) %>%
  as.data.frame(.),
data = df_wide %>%
  # select(-yday, -year, -alb)%>%
  select(-date, -year, -yday) %>%
  as.data.frame(.), # %>%
# tail(1000),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 3))
)


#
#
#
# #do we really have to list out all the vars in the formula?
# # or can we just ditch the ones we don't care about
# pred_sdl = mnimput(df_wide%>%
#                      select(-year)%>%
#                      as.data.frame(.), data = df_wide%>%
#                      #select(-yday, -year, -alb)%>%
#                      select(-year)%>%
#                      as.data.frame(.),#%>%
#                    #tail(1000),
#                    maxit =50, ts = TRUE,
#                    log = TRUE, log.offset = 0.01)


# t2= predict(tst)
sdl_infl <- predict(pred_sdl)
sdl_infl$date <- df_wide %>% # tail(1000)%>%
  dplyr::pull(date)

sdl_infl <- sdl_infl %>%
  select(date, mar, sdl) %>%
  tidyr::pivot_longer(cols = mar:sdl, names_to = "local_site", values_to = "discharge")

# quick plot of how the infilling looks
if (plot_in_code) {
  ggplot(sdl_infl %>% filter(local_site == "sdl" &
    date >= "1998-10-01"), aes(x = date, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(color = "red") +
    geom_point(data = infilled_df %>% filter(local_site == "sdl"), aes(x = date, y = discharge), color = "blue") +
    # geom_line(color="blue")+
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}

# bring the saddle predictions in
infilled_df <- infilled_df %>%
  # only infill starting in winter of 1998, spring 1999 is the earliest year
  # with any saddle data
  left_join(sdl_infl %>% filter(local_site == "sdl" & date >= "1998-10-01") %>%
    rename(sdl_disch_inf = discharge))

if (plot_in_code) {
  # looks ok-ish, though
  # might consider xts infilling for gaps <2d?
  ggplot(infilled_df %>% filter(local_site == "sdl" &
    date >= "1998-10-01"), aes(x = yday, y = discharge)) +
    geom_point(
      data =
        infilled_df %>% filter(local_site == "sdl" &
          date >= "1998-10-01" & is.na(discharge)),
      aes(y = sdl_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)

  # this looks ok, we could potentially
  # add in some 1999 off season if wanted one more year
  ggplot(infilled_df %>% filter(local_site == "sdl"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "sdl" & is.na(discharge)),
      aes(y = sdl_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}



# repeat infilling mart off sdl
# use mnimput to predict the missing mar discharge
# does this look ok
df_wide <- infilled_df %>%
  filter(yday > min(flow_season$start[flow_season$local_site %in% c("sdl", "mar")]) &
    yday < max(flow_season$end[flow_season$local_site %in% c("sdl", "mar")])) %>%
  filter(local_site %in% c("mar", "sdl")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
  arrange(date)

# add in the mar for all years
for (yr in unique(df_wide$year)) {
  newdat <- df %>%
    filter(local_site == "mar" & lubridate::year(date) == yr) %>%
    mutate(yday = lubridate::yday(date))
  if (yr == 2013) {
    newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
  }
  # skip as predicators any years that are all na
  if (!all(is.na(newdat$discharge))) {
    newdat[[paste0("mar_", yr)]] <- newdat[["discharge"]]
    df_wide <- df_wide %>%
      left_join(., newdat %>%
        select(
          -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
          -notes
        ))
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$year)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$year == yr] <- NA
  }
}

df_wide <- df_wide %>%
  # remove flood
  mutate(mar = ifelse(year == 2013 & yday %in% c(253:259), NA, mar)) %>%
  arrange(date)
# make predictions

pred_mar <- mnimput(df_wide %>%
  select(-date, -year, -yday) %>%
  as.data.frame(.),
data = df_wide %>%
  # select(-yday, -year, -alb)%>%
  select(-date, -year, -yday) %>%
  as.data.frame(.), # %>%
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 3))
)


mar_infl <- predict(pred_mar)
mar_infl$date <- df_wide %>% # tail(1000)%>%
  dplyr::pull(date)

mar_infl <- mar_infl %>%
  select(date, mar, sdl) %>%
  tidyr::pivot_longer(cols = mar:sdl, names_to = "local_site", values_to = "discharge")

if (plot_in_code) {
  # vis check
  ggplot(mar_infl %>% filter(local_site == "mar" & date < "2019-10-07"), aes(x = date, y = discharge)) +
    geom_line(data = infilled_df %>% filter(local_site == "mar"), aes(x = date, y = discharge), color = "red") +
    geom_point(data = infilled_df %>% filter(local_site == "mar"), aes(x = date, y = discharge), color = "red") +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}


# end infilling mart off sdl
# add to df
infilled_df <- infilled_df %>%
  left_join(mar_infl %>% filter(local_site == "mar" & date < "2019-10-07") %>%
    rename(mar_disch_inf = discharge))
if (plot_in_code) {
  # looks ok except 1987 when started too late?
  ggplot(infilled_df %>% filter(local_site == "mar" & date < "2019-10-07"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "mar" &
          date < "2019-10-07" & is.na(discharge)),
      aes(y = mar_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)


  ggplot(infilled_df %>% filter(local_site == "mar" & date < "2019-10-07" & wyear == 1987), aes(x = yday, y = discharge)) +
    geom_point(
      data =
        infilled_df %>% filter(local_site == "mar" &
          date < "2019-10-07" & is.na(discharge) & wyear == 1987),
      aes(y = mar_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}



### gl4 infilling###############################################
# infill gl4 off self in other wyrs (correlation with other sites basically
# goes to 0 in winter months, which is what we are missing)

df_wide <- infilled_df %>%
  filter(local_site %in% c("gl4")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>% # ,
  arrange(date) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1))

# add in the gl4 for all wyears
for (yr in unique(df_wide$wyear)) {
  if (yr == 2021) {
    next
  } else { # 1 day only
    if (yr == 2015) {
      next
    } else { # sensor appears bad
      newdat <- df %>%
        mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
        # month = lubridate::month(date)) %>%
        mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
        filter(local_site == "gl4" & wyear == yr)
      if (yr == 2013) {
        newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
      }

      # skip as predicators any years that are all na
      if (!all(is.na(newdat$discharge))) {
        newdat[[paste0("gl4_", yr)]] <- newdat[["discharge"]]
        df_wide <- df_wide %>%
          left_join(., newdat %>%
            select(
              -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
              -notes, -year, -wyear
            ))
      }
    }
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$wyear)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$wyear == yr] <- NA
  }
}

# do not use correlations with wyear 2015 in the predictions
# as these data are bad, similarly
# we expect the correlations with the flood data across years not to be relevant
df_wide <- df_wide %>%
  arrange(date) %>%
  mutate(gl4 = ifelse(wyear == 2015, NA, gl4)) %>%
  # sce should do this for mar, sad too
  mutate(gl4 = ifelse(wyear == 2013 & yday %in% c(253:259), NA, gl4))

# make predictions
pred_gl4 <- mnimput(df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
data = df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 4))
)

gl4_infl <- predict(pred_gl4)
gl4_infl$date <- df_wide %>%
  dplyr::pull(date)

gl4_infl <- gl4_infl %>%
  select(date, gl4) %>%
  tidyr::pivot_longer(cols = gl4, names_to = "local_site", values_to = "discharge")

if (plot_in_code) {
  # vis check
  ggplot(gl4_infl %>% filter(local_site == "gl4"), aes(x = date, y = discharge)) +
    geom_line(color = "blue") +
    geom_line(
      data = infilled_df %>% filter(local_site == "gl4"),
      aes(x = date, y = discharge), color = "red"
    ) +
    facet_wrap(~local_site, scales = "free_y", ncol = 1)
}

# join to the rest of the data
infilled_df <- infilled_df %>%
  left_join(gl4_infl %>% filter(local_site == "gl4" & date < "2020-10-01") %>%
    rename(gl4_disch_inf = discharge)) %>%
  # remove wy 2015 seems broken sensor of some sort and 2018 hardly any data
  mutate(
    gl4_disch_inf = ifelse(wyear %in% c(2015, 2018), NA, gl4_disch_inf),
    discharge = ifelse(wyear == 2015 & local_site == "gl4", NA, discharge)
  )

if (plot_in_code) {
  # looks okish?
  ggplot(infilled_df %>% filter(local_site == "gl4" & date < "2020-10-01"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "gl4" &
          date < "2020-10-01" & is.na(discharge)),
      aes(y = gl4_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y", ncol = 4)
}
# repeat for albion##############################################################
### gl4 infilling###############################################
# infill alb off self in other wyrs (correlation with other sites basically
# goes to 0 in winter months, which is what we are missing)

df_wide <- infilled_df %>%
  filter(local_site %in% c("alb")) %>%
  select(local_site, date, discharge) %>%
  tidyr::pivot_wider(names_from = local_site, values_from = discharge, values_fill = NA) %>%
  mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>% # ,
  arrange(date) %>%
  mutate(wyear = ifelse(yday < 274, year, year + 1))

# add in the alb for all wyears
for (yr in unique(df_wide$wyear)) {
  if (yr == 2021) {
    next
  } else { # 1 day only
    newdat <- df %>%
      mutate(yday = lubridate::yday(date), year = lubridate::year(date)) %>%
      mutate(wyear = ifelse(yday < 274, year, year + 1)) %>%
      filter(local_site == "alb" & wyear == yr)
    if (yr == 2013) {
      newdat$discharge[newdat$yday %in% c(253:259)] <- NA # remove 2013 flood as seasonal predictor
    }

    # skip as predicators any years that are all na
    if (!all(is.na(newdat$discharge))) {
      newdat[[paste0("alb_", yr)]] <- newdat[["discharge"]]
      df_wide <- df_wide %>%
        left_join(., newdat %>%
          select(
            -file, -date, -LTER_site, -local_site, -date, -discharge, -temperature,
            -notes, -year, -wyear
          ))
    }
  }
}

# remove the self-self predictions
for (yr in unique(df_wide$wyear)) {
  mycol <- names(df_wide)[grepl(yr, names(df_wide))]
  if (length(mycol) > 0) {
    df_wide[[mycol]][df_wide$wyear == yr] <- NA
  }
}

# do not use correlations with wyear 2015 in the predictions
# as these data are bad, similarly
# we expect the correlations with the flood data across years not to be relevant
df_wide <- df_wide %>%
  arrange(date) %>%
  mutate(alb = ifelse(wyear == 2013 & yday %in% c(253:259), NA, alb))

# make predictions
pred_alb <- mnimput(df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
data = df_wide %>%
  select(-date, -year, -yday, -wyear) %>%
  as.data.frame(.),
maxit = 50, ts = TRUE,
log = TRUE, log.offset = 0.01,
sp.control =
  list(df = rep(sp_control, ncol(df_wide) - 4))
)

alb_infl <- predict(pred_alb)
alb_infl$date <- df_wide %>%
  dplyr::pull(date)

alb_infl <- alb_infl %>%
  select(date, alb) %>%
  tidyr::pivot_longer(cols = alb, names_to = "local_site", values_to = "discharge")

# join to the rest of the data
infilled_df <- infilled_df %>%
  left_join(alb_infl %>% filter(local_site == "alb" &
    date <= max(df$date[!is.na(df$discharge) & df$local_site == "alb"])) %>%
    rename(alb_disch_inf = discharge))

if (plot_in_code) {
  # looks ok-ish except 1999 is clearly
  # all just super low (some management?)
  ggplot(infilled_df %>% filter(local_site == "alb" & date < "2020-10-01"), aes(x = yday, y = discharge)) +
    # geom_line(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='red')+
    # geom_point(data = infilled_df%>%filter(local_site=='sdl'), aes(x=date, y=discharge), color='blue')+
    geom_point(
      data =
        infilled_df %>% filter(local_site == "alb" &
          date < "2020-10-01" & is.na(discharge)),
      aes(y = alb_disch_inf), color = "red"
    ) +
    geom_point(color = "blue") +
    facet_wrap(~wyear, scales = "free_y")
}
# end infilling


bu <- infilled_df

##### join infilling to main data and plot all sites
# plot all
infilled_df <- infilled_df %>%
  mutate(is_infilled = ifelse(is.na(discharge), TRUE, FALSE)) %>%
  mutate(is_infilled = ifelse(yday < start | yday > end, TRUE, is_infilled)) %>%
  # gaps of up to 3 d estimate by linear interpolation
  group_by(local_site) %>%
  arrange(local_site, date) %>%
  mutate(discharge = zoo::na.approx(discharge, maxgap = 3, method = "linear")) %>%
  ungroup() %>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "mar",
    mar_disch_inf, discharge
  )) %>%
  # wyear<2020, mar_disch_inf, discharge))%>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "sdl",
    sdl_disch_inf, discharge
  )) %>% # &
  ## wyear>1998, sad_disch_inf, discharge))
  mutate(discharge = ifelse(is.na(discharge) & local_site == "gl4" & !wyear %in% c(2015, 2018),
    gl4_disch_inf, discharge
  )) %>%
  mutate(discharge = ifelse(is.na(discharge) & local_site == "alb",
    alb_disch_inf, discharge
  ))


infilled_df <- infilled_df %>%
  # gaps of up to 3 d in mart/sdl againn estimate by linear interpolation
  group_by(local_site) %>%
  arrange(local_site, date) %>%
  mutate(discharge = zoo::na.approx(discharge, maxgap = 3, method = "linear")) %>%
  ungroup()


# handful of very slightly negative values from
# the back transform

infilled_df <- infilled_df %>%
  mutate(discharge = ifelse((!is.na(discharge) & discharge < 0), 0, discharge))


write.csv(infilled_df %>%
  select(local_site, date, discharge, yday, year, wyear, is_infilled),
paste0("infilled_data/spctl_", sp_control, ".csv"),
row.names = FALSE
)

#infilled_df = read.csv("discharge/infilled_data/spctl_7.csv")

#calc cum curves


(ggplot(
  infilled_df,
  aes(x = date, y = discharge, color = is_infilled)
) +
  geom_point(size = 0.01) +
  facet_wrap(~local_site, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "years") +
  theme(axis.text.x = element_text(angle = 90))) %>%
  ggsave(
    file = paste0("plots/overtime_spctl_", sp_control, ".jpg"),
    width = 8, height = 4
  )

(ggplot(
  infilled_df %>% filter(is_infilled == FALSE),
  aes(x = yday, y = discharge, color = factor(year))
) +
  geom_line(size = 0.01) +
  facet_wrap(~local_site, scales = "free_y", ncol = 1) +
  # scale_x_date(date_breaks = "years")+
  theme(axis.text.x = element_text(angle = 90))) %>%
  ggsave(
    file = paste0("plots/raw_by_year.jpg"),
    width = 8, height = 4
  )




for (loc in unique(infilled_df$local_site)) {
  (ggplot(
    infilled_df %>% filter(local_site == loc) %>%
      mutate(yday = ifelse(wyear != year, yday - 365, yday)),
    aes(x = yday, y = discharge, color = is_infilled)
  ) +
    geom_point(size = 0.01) +
    facet_wrap(~wyear, ncol = 2) +
    # scale_x_date(date_breaks = "years")+
    theme(axis.text.x = element_text(angle = 90))) %>%
    ggsave(
      file = paste0("plots/", loc, "byyear_spctl_", sp_control, ".jpg"),
      height = 20, width = 5
    )
}

# sce stop loop here
# consider weighting by taking - over all yrs, mean proportion of total
# discharge on each doy, then taking the weighted sum of days missing for wts,
# if 0 data missing wts = 1, if all data missing wts = 0, if 50% of discharge
# infilled wts = 0.5
# if 25% of discharge missing wt = 0.75


# also consider comparing all 4 versions to see if they look any different?


wts <- infilled_df %>%
  filter(is_infilled == FALSE) %>%
  group_by(local_site, yday) %>%
  summarize(mean_daily = mean(discharge, na.rm = TRUE)) %>%
  left_join(infilled_df %>%
    filter(is_infilled == FALSE) %>%
    group_by(local_site, yday) %>%
    summarize(mean_daily = mean(discharge, na.rm = TRUE)) %>%
    filter(yday != 366) %>%
    group_by(local_site) %>%
    summarize(tot = sum(mean_daily))) %>%
  mutate(prop = mean_daily / tot) %>%
  full_join(., infilled_df)


wts <- wts %>%
  # always remove yday 366 so all yrs have same number of days
  filter(yday < 366) %>%
  # mutate(prop_this_year = prop * discharge)%>%
  group_by(local_site, wyear, is_infilled) %>%
  summarise(tot_disch = sum(prop, na.rm = TRUE)) %>%
  # summarise(tot_disch = sum(prop_this_year, na.rm=TRUE))%>%
  pivot_wider(names_from = is_infilled, values_from = tot_disch) %>%
  mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(prop_infilled = `TRUE` / (`FALSE` + `TRUE`))



# find max date by wyear
peak_flow <- infilled_df %>%
  group_by(wyear, local_site) %>%
  mutate(disch_ma_2 = slider::slide_dbl(discharge, mean,
    .before = 2,
    .after = 2
  )) %>%
  summarize(max_disch = max(disch_ma_2, na.rm = TRUE)) %>%
  full_join(infilled_df) %>%
  mutate(disch_ma_2 = slider::slide_dbl(discharge, mean,
    .before = 2,
    .after = 2
  )) %>%
  rowwise() %>%
  filter(max_disch == disch_ma_2) %>%
  filter(wyear < 2021)

# need to pull out some weird albion winter floods
# here, prob the 2013 floods too, and the last year
# but there's not really a super strong movement of the peak date
# more that the q50 is getting earlier but not the peak
ggplot(peak_flow, aes(x = wyear, y = yday)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site) +
  ggtitle("timing of peak flow")

# do q50s for these
tots <- infilled_df %>%
  # always remove yday 366 so all yrs have same number of days
  filter(yday < 366) %>%
  group_by(local_site, wyear) %>%
  summarise(
    tot_disch = sum(discharge, na.rm = TRUE),
    tot_day = sum(!is.na(discharge))
  ) %>%
  # use only complete years
  mutate(tot_disch = ifelse(tot_day != 365, NA, tot_disch)) %>%
  mutate(
    disch_50 = tot_disch / 2,
    disch_20 = tot_disch / 5,
    disch_80 = 4 * (tot_disch / 5)
  ) %>%
  left_join(., wts %>% select(local_site, wyear, prop_infilled))





#
# tots = infilled_df %>%
#   #always remove yday 366 so all yrs have same number of days
#   filter(yday<366)%>%
#   group_by(local_site, wyear)%>%
#   summarise(tot_disch = sum(discharge, na.rm=TRUE),
#             tot_day = sum(!is.na(discharge)))%>%
#   #use only complete years
#   mutate(tot_disch = ifelse(tot_day!=365, NA, tot_disch))%>%
#   mutate(disch_50 = tot_disch/2,
#          disch_20 = tot_disch/5,
#          disch_80 = 4*(tot_disch/5))









avg_tots <- tots %>%
  group_by(local_site) %>%
  summarize(tot_disch = mean(tot_disch, na.rm = TRUE)) %>%
  mutate(
    disch_50_all = tot_disch / 2,
    disch_20_all = tot_disch / 5,
    disch_80_all = 4 * (tot_disch / 5)
  )

#        nday = )%>%
# cut out incomplete years by site
# mutate(tot_disch = ifelse((local_site=='mar'&(wyear>2019|wyear<1983)), NA, tot_disch))%>%
# mutate(tot_disch = ifelse((local_site=='sdl'&(wyear>2020|wyear<2000)), NA, tot_disch))%>%
# mutate(tot_disch = ifelse((wyear>2020), NA, tot_disch))%>%
# mutate(tot_disch = ifelse((wyear%in%c(2015, 2018)&local_site=='gl4'), NA, tot_disch))%>%
# mutate(tot_disch = ifelse((wyear%in%c(1981)&local_site=='alb'), NA, tot_disch))%>%
# mutate(disch_50 = tot_disch/2,
#        disch_20 = tot_disch/5,
#        disch_80 = 4*(tot_disch/5))

# plot
q50 <- infilled_df %>%
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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

# sce need to remake these so the q20 can be negative if it's in the fall
# for alb etc.
q20_mean <- infilled_df %>%
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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
  # select(-tot_disch)%>%
  # left_join(., tots)%>%
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
  mutate(duration = q80 - q20) %>%
  mutate(duration_mean = q80_mean - q20_mean)



# martinelli does seem like it's discharging a bit earlier over time.
# but over just the past 20 yrs
# discharge is actually a bit later in both earlier over time.
ggplot(q50, aes(x = wyear, y = q50, color = prop_infilled)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site) +
  ggtitle("q50 over time")

# alb and gl4 to melt out earlier (spring)
# but note we need to pee


ggplot(q50, aes(x = wyear, y = q50_mean, color = prop_infilled)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site) +
  ggtitle("q50 (date hitting all-time mean 50% discharge) over time")



# strong relationship as expected for snowmelt-dominated systems
# the median discharge date is later in high total discharge years
ggplot(q50, aes(x = tot_disch, y = q50)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site, scales = "free_x")


# but they hit the 50% of all-time average discharge earlier in high discharge
# yrs bc the whole discharge is elevated
ggplot(q50, aes(x = tot_disch, y = q50_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site, scales = "free_x")

# which I think means in high water years vs low water years, the most
# pronounced differences are in the 2nd half of the season? ie
# water is pretty invariant in the spring but more variable in the summer drying
# out period

# duration stuff starts here, not quite right
vars <- infilled_df %>%
  filter(yday != 366) %>%
  mutate(month = lubridate::month(date)) %>%
  group_by(local_site, year, month) %>%
  summarize(tot = sum(discharge), nday = dplyr::n()) %>%
  filter(nday >= 28)

# for gl4 and albion
# years with more tot discharge, the 60% of middle flows
# are spread out over a shorter period of time (q20 - q80)
# this is not the case for mart or sdl though
ggplot(q50, aes(x = tot_disch, y = duration)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site, scales = "free_x")

# if you look at the seasonality of the average 20-80 pct flows
# for mart and saddle, the high discharge years are actually the shortest
ggplot(q50, aes(x = tot_disch, y = duration_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site, scales = "free_x")


# show the relationship, earlier years, longer duration
# in general..
ggplot(q50, aes(x = q50, y = duration_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~local_site)

# we really need to add peak swe to this
# but
mod <- lm(duration ~ tot_disch + q50,
  data =
    q50 %>% filter(local_site == "mar")
)


#infilling of gl4, alb start HERE

