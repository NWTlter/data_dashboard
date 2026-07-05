# plot sdl spp comp by growth form
# and NPP
# SCE 10 Jan 2022

# Credit for basic workflow for downloading/top subsetting/cleaning from
# M Oldfather script in nwt_8-renewal/get_saddle_sp_CNM.R

# -- SETUP -----
# clean up enviro, read in needed libraries
rm(list=ls())
#library(readxl)
library(tidyverse)
#library(vegan)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utility_functions/utility_functions_all.R")

#TODO - need a few updates to the published taxon list then can download
# from EDI


# download data -----------------------------------------------------------

#saddle grid spp comp
sdlcomp <- getTabular(93)

#saddle grid ANPP
sdlprod <- getTabular(16)

#get growth forms
#sdlgf <- getTabular(91)

#need to get this onto EDI with the growth forms, but
# use version emailed from Jane-> Jared for now hopefully correct
sdlgf <-read.csv('data/pspecies.mw.data_JGS.csv', na.strings = na_vals)


# Saddle comp by growth form ----------------------------------------------

#subset to usda only not duplicated
gf_join = sdlgf %>%
  select(USDA_code, USDA_family, category,group, growth_habit) %>%
  distinct() %>%
  filter(!duplicated(USDA_code)&!duplicated(USDA_code), fromLast = TRUE)

sdlcomp <-sdlcomp %>%
  left_join(., gf_join)
#only one missing is Noccaea montana and a bunch of carexes
unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)|sdlcomp$growth_habit==""])
unique(sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])

#could more efficiently do with case_when but ordered mutate works ok for now
sdlcomp = sdlcomp %>%
  mutate(growth_habit = ifelse((is.na(growth_habit)&grepl('Rock|Bare|Hole|soil', USDA_name)), 'soil', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Scat|scat', USDA_name)), 'scat', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Carex', USDA_name)), 'graminoid', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Carex', USDA_name)), 'graminoid', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Mert|Nocc|Erys|Ante|compos|Pack', USDA_name)), 'forb', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Litter', USDA_name)), 'litter', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('Unknown', USDA_name)), 'unknown', growth_habit),
         growth_habit = ifelse((is.na(growth_habit)&grepl('marker', USDA_name)), 'marker', growth_habit))

# unique (sdlcomp$USDA_name[is.na(sdlcomp$growth_habit)])
# View (sdlcomp %>% filter())

# change year = 1996 to 1995 for plot 37 (sampled one year late for a couple of plots)
sdlcomp[sdlcomp$year == 1996, "year"] <- 1995

# subset to only top hit (or bottom hit if no top hit present), by point
sdl_top<-
  sdlcomp %>% 
  filter(hit_type == "bottom" | hit_type == "top") %>% 
  arrange(year, plot, x, y, desc(hit_type)) %>% 
  group_by(year, plot, x,y) %>% 
  slice(1) %>%
  ungroup()%>%
  #add subshrubs to shrubs (Salix nivalis)
  mutate(growth_habit = ifelse(grepl('shrub', growth_habit), 'shrub', growth_habit))

#sdl_top

# aggregate data_file to calculate the number of hits
number_hits <- sdl_top %>% 
  #add subshrubs to shrubs (Salix nivalis)
  group_by(year, plot, growth_habit) %>% 
  summarise(hits = dplyr::n()) 


# Fill in zeros for species/non-species not hit in a certain plot/year
data_with_zeros <- expand_grid(year = unique(number_hits$year), plot = unique(number_hits$plot), growth_habit = unique(number_hits$growth_habit))
data_with_zeros

# make abundance data frame
number_hits <- left_join(data_with_zeros, number_hits)
number_hits[is.na(number_hits$hits), "hits"] <- 0


# aggregate w the zeros
number_hits_by_year <- number_hits%>% 
  #add subshrubs to shrubs (Salix nivalis)
  group_by(year, growth_habit) %>% 
  summarise(mean_cover = mean(hits)) 

ggplot (number_hits_by_year,
        aes(x=year, y=mean_cover, group = growth_habit, color = growth_habit))+
  geom_line()+
  facet_wrap(~growth_habit)


anom_growth_form <- number_hits_by_year%>%
  group_by(growth_habit) %>%
  mutate(
    avgcov = mean(mean_cover, na.rm = TRUE),
    anom_cov = mean_cover - avgcov,
    posneg = ifelse(anom_cov > 0, "pos", "neg") %>%
      factor(c("pos", "neg")),
    growth_habit = ifelse(growth_habit == "nonvascular", 'moss', growth_habit),
    growth_habit = ifelse(growth_habit == "lichenous", 'lichen', growth_habit)
  ) %>%
  filter(growth_habit %in% c('forb', 'graminoid',
                             'shrub', 'moss', 'lichen')) %>%
  mutate(growth_habit = factor(growth_habit,
                               levels = c('forb', 'graminoid',
                               'shrub', 'moss', 'lichen')))

write.csv(anom_growth_form, 'data_derived/sdl_veggf_anom.csv', row.names = FALSE)
           
# plot
g1 <- ggplot(anom_growth_form , aes(x = year, y = anom_cov)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  labs(y = "Percent cover \n difference from long-term mean", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(more %<->% less))) +
  theme_hc() +
  facet_wrap(~growth_habit, scales = 'free_y')+
  theme(legend.position = "none")

ggsave(g1,
       file = "figs/sdl_gf_anom.png",
       scale = 0.8, width = 10, height = 6
)



# NPP ---------------------------------------------------------------------

npp_by_year <- sdlprod%>% 
  #aggregat from when n=2 subsamples
  group_by(year, grid_pt) %>% 
  summarise(NPP = mean(NPP)) %>%
  #group_by(veg_class)
  ungroup() %>%
  group_by(year) %>%
  summarise(NPP = mean(NPP)) %>%
  ungroup()%>%
  mutate(
    avgNPP = mean(NPP, na.rm = TRUE),
    anom_NPP = (NPP*100/ avgNPP)-100,
    posneg = ifelse(anom_NPP > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

# plot NPP
g1 <- ggplot(npp_by_year , aes(x = year, y = anom_NPP)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  labs(y = "Above ground biomass \n difference from long-term mean", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = I, breaks = NULL, name = expression(more %<->% less))) +
  theme_hc() +
  theme(legend.position = "none")

ggsave(g1,
       file = "figs/sdl_npp_anom.png",
       scale = 0.8, width = 10, height = 6
)

