# Script to plot c02 trends for data dashbord
# SCE updated 12 Aug 2025

library(tidyverse)
library(ggplot2)
library(padr)
library(scales)
library (ggthemes)

# understand metadata -----------------------------------------------------
# https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/txt/co2_mlo_surface-flask_1_ccgg_event.txt

# citation and readme here
# https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/README_co2_surface-flask_ccgg.html

nwt <- read.table("https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/txt/co2_nwr_surface-flask_1_ccgg_month.txt") %>%
  mutate(date = as.Date(paste0("01/", V3, "/", V2), format = "%d/%m/%Y")) %>%
  rename(CO2 = V4, year = V2, month = V3)

# use padr to pad dates
nwt <- nwt %>%
  padr::pad()

nwt_event <- read.table("https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/txt/co2_nwr_surface-flask_1_ccgg_event.txt")

names(nwt_event) <- nwt_event[1, ]
nwt_event <- nwt_event[-1, ]


# notes on qc flags
# other than . in first column is problems with analysis
# . in second doesn't meet some criterion for analysis (ok for this)
# . p in third column is prelim
# so we'll make sure there's a . starting the 2nc column

nwt_event <- nwt_event %>%
  filter(grepl("^\\.", qcflag)) %>%
  mutate(
    datetime = lubridate::ymd_hms(datetime),
    value = as.numeric(value),
    # extract just the ymd component of datetime
    date = as.Date(datetime)
  ) %>%
  rename(CO2 = value)


# plotting ----------------------------------------------------------------


g1 <- ggplot(
  nwt,
  aes(x = date, y = CO2)
) +
  geom_point(
    data = nwt_event, aes(x = date, y = CO2),
    color = "red", size = 0.001, alpha = 0.2
  ) +
  geom_line(size = 0.1, alpha = 1) +
  ylab(bquote(CO[2] ~ (mu * mol ~ mol^{
    -1
  }))) +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  scale_x_date(
    breaks = seq(as.Date("1960-01-01"), as.Date("2030-01-01"), by = "10 years"),
    date_labels = "%Y",
    limits = c(as.Date("1960-01-01"), Sys.Date())
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  xlab("Year") +
  geom_smooth(
    data = nwt_event, aes(x = date, y = CO2), method = "loess", # span = 0.2,
    se = TRUE, color = "transparent", linetype = "solid",
    fill = "blue"
  ) +
  theme_hc()


ggsave(g1,
  file = "noaa_gml/figures/noaa_co2.png",
  scale = 0.5, width = 10, height = 6
)
