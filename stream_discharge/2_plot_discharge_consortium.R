# plot discharge stats off infilled data
# SCE 10 Jan 2023


# -- SETUP -----
# clean up enviro, read in needed libraries
rm(list=ls())
library(tidyverse)
library (lemon)
library (ggthemes)
library (lubridate)
library (ggplotFL)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read gap-filled data, created in script q_50.r
# takes a while to run so not redone here

#span control had little effect on infilling, so let's pick 7
infilled_df <-read_csv("infilled_data/spctl_7.csv")


#calc cum curves
# (ggplot(
#   infilled_df,
#   aes(x = date, y = discharge, color = is_infilled)
# ) +
#     geom_point(size = 0.01) +
#     facet_wrap(~local_site, scales = "free_y", ncol = 1) +
#     scale_x_date(date_breaks = "years") +
#     theme(axis.text.x = element_text(angle = 90))) %>%
#   ggsave(
#     file = paste0("plots/overtime_spctl_", sp_control, ".jpg"),
#     width = 8, height = 4
#   )
# 
# (ggplot(
#   infilled_df %>% filter(is_infilled == FALSE),
#   aes(x = yday, y = discharge, color = factor(year))
# ) +
#     geom_line(size = 0.01) +
#     facet_wrap(~local_site, scales = "free_y", ncol = 1) +
#     # scale_x_date(date_breaks = "years")+
#     theme(axis.text.x = element_text(angle = 90))) %>%
#   ggsave(
#     file = paste0("plots/raw_by_year.jpg"),
#     width = 8, height = 4
#   )

# for (loc in unique(infilled_df$local_site)) {
#   (ggplot(
#     infilled_df %>% filter(local_site == loc) %>%
#       mutate(yday = ifelse(wyear != year, yday - 365, yday)),
#     aes(x = yday, y = discharge, color = is_infilled)
#   ) +
#     geom_point(size = 0.01) +
#     facet_wrap(~wyear, ncol = 2) +
#     # scale_x_date(date_breaks = "years")+
#     theme(axis.text.x = element_text(angle = 90))) %>%
#     ggsave(
#       file = paste0("plots/", loc, "byyear_spctl_", sp_control, ".jpg"),
#       height = 20, width = 5
#     )
# }


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

# are the wyears complete
cts = infilled_df%>%
  group_by(local_site, wyear) %>%
  summarize(n = dplyr::n())

# yes except albion 1981 an da few
# days of sld in the last year after it was already 0


# make all into anomalies, by month or not
anom_disch <- infilled_df%>%
  filter(wyear >1981&wyear<2021) %>%
group_by(local_site, wyear) %>%
  summarise(
    total_disch = sum(discharge),
    ct = dplyr::n())%>%
  ungroup()%>%
  group_by(local_site) %>%
  mutate(avgdisch = mean(total_disch, na.rm = TRUE),
    anom_disch = total_disch *100/avgdisch -100,
    posneg = ifelse(anom_disch > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))) %>%
  ungroup()


facet_names <-  c(alb = "Albion",
                  gl4 = "Green Lake 4",
                  mar = "Martinelli",
                  sdl = "Saddle")

g1 <- ggplot(anom_disch , aes(x = wyear, y = anom_disch)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#034e7b", "#99000d")) +
  labs(y = "Stream discharge \n difference from long-term mean (%)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(more~streamflow %<->% less~streamflow))) +
  theme_hc() +
  facet_wrap(~local_site, scales = 'free_y',
             labeller = as_labeller(facet_names))+
  theme(legend.position = "none")


ggsave(g1,
       file = "plots/streamflow_anom.png",
       scale = 0.6, width = 10, height = 6
)


# do just gl4 by month
# make all into anomalies, by month or not
# should be at least 30ish days
anom_disch_month_gl4 <- infilled_df%>%
  filter(local_site == 'gl4'&wyear<2021) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(5:9)) %>%
  group_by(wyear, month) %>%
  summarise(
    total_disch = sum(discharge))%>%
  ungroup()%>%
  group_by(month) %>%
  mutate(avgdisch = mean(total_disch, na.rm = TRUE),
         anom_disch = total_disch *100/avgdisch -100,
         posneg = ifelse(anom_disch > 0, "pos", "neg") %>%
           factor(c("pos", "neg"))) %>%
  ungroup() %>%
  mutate(month_name = 
           month(month, label=TRUE))


#seems hard to view this way
g2 <- ggplot(anom_disch_month_gl4, aes(x = wyear, y = anom_disch)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("#034e7b", "#99000d")) +
  labs(y = "Stream discharge \n difference from long-term mean (%)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(more~streamflow %<->% less~streamflow))) +
  theme_hc() +
  facet_wrap(~month_name)+
  theme(legend.position = "none")

ggsave(g2,
       file = "plots/streamflow_anom_gl4_month.png",
       scale = 0.7, width = 10, height = 6
)

summer_disch_last_yr <- infilled_df%>%
  filter(local_site == 'gl4'&wyear<2021) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(5:9)) %>%
  mutate(yday = lubridate::yday(date)) %>%
  filter(wyear == 2020) %>%
  left_join(infilled_df%>%
              filter(local_site == 'gl4'&wyear<2021) %>%
              mutate(yday = lubridate::yday(date),
                     month = lubridate::month(date)) %>%
              filter(month %in% c(5:9)) %>%
              filter(wyear != 2020) %>%
              select(yday, discharge) %>%
              rename(oth = discharge))

g3 <- ggplot(summer_disch_last_yr %>%
               ungroup(), aes(x = date)) +
  geom_flquantiles(aes(y = oth),
                   probs = c(0.10, 0.90), fill = "red", alpha = 0.25,
                   na.rm = TRUE
  ) +
  geom_line(aes(y = discharge)) +
  ylab("Streamflow (cubic meters)") +
  theme_hc()+
  ggtitle('Green Lake 4 2020 \n summer streamflow seasonality \n vs. historical patterns')

ggsave(g3,
       file = "plots/gl4_summer_disch_2020.png",
       scale = 0.6, width = 8, height = 6
)




