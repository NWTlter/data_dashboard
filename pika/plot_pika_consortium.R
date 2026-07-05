# Created by SCE 1/12/2022
# Based on Cray equivalent script
# that has more details but different aesthetics in figs
# script-pika-juvie-capture-ratio-vs-GDD.R
##########################

# -- SETUP ------
# libraries
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# TODO - set up to run directly from EDI and update

# -- READ DATA ------
# current pika demog + older cleared of recapture, combined
d <- read.csv("C:/Users/Sarah/Documents/git/nwt8-renewal/final_figs/pika/inputs/pika-demog-WK-1981-2020.csv") # WK pikas only

# prep pika data

dt <- strptime(d$Date, "%m/%d/%Y")

d <- mutate(d, day = dt$mday, month = dt$mon + 1, year = dt$year + 1900, doy = dt$yday, count = 1)

by.yr <- summarize(group_by(d, year),
  juvies = sum(count[Stage == "J"]),
  adults = sum(count[Stage == "A"]),
  j.mu = mean(Wt.g[Stage == "J"], na.rm = T),
  j.sd = sd(Wt.g[Stage == "J"], na.rm = T),
  j.n = sum(!is.na(Wt.g[Stage == "J"]), na.rm = T),
  a.n = sum(!is.na(Wt.g[Stage == "A"]), na.rm = T),
  capdate.mu = mean(doy, na.rm = T)
)


# the high juvie-to-adult ratio in 2018 was an outlier that probably should be
# censored because it was based on just 2 adults (!), the lowest ever captured
p1 <- ggplot() +
  geom_point(
    data = by.yr %>%
      filter(year != 2018),
    aes(year, y = juvies / adults),
    size = 3
  )

# make as anomalies
by.yr.anom <- by.yr %>%
  mutate(j.a = juvies / adults) %>%
  mutate(
    meanja = mean(j.a, na.rm = TRUE),
    ja_anom = j.a * 100 / meanja - 100,
    posneg = ifelse(ja_anom > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )


g1 <- ggplot(by.yr.anom, aes(x = year, y = ja_anom)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#034e7b", "#99000d")) +
  labs(y = "Pika fecundity anomaly (%) \n ratio of juveniles:adults relative to long-term mean", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(
    trans = I,
    breaks = NULL, name = expression(high ~ reproduction %<->% low ~ reproduction)
  )) +
  theme_hc() +
  theme(legend.position = "none")

ggsave(g1,
  file = "plots/pika_ja_anom.png",
  scale = 0.7, width = 8, height = 6
)
