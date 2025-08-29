# NWT LTER Saddle Composition Data Analysis
# 20210601

library(tidyverse)
library(ggplot2)
library(lme4)
library(MASS)
library(vegan)
library(picante)
library(codyn)
library(cowplot)
library(picante)
library(ape)
library(ggforce)
library(lubridate)
library(lmerTest)
library(emmeans)
library(viridisLite)
library(geomtextpath)
library(ggrepel)
library(directlabels)
library(here)

# code to get most recent saddle data (not currently working)
#source("utility_functions/utility_functions_all.R")
#data_file <-getTabular(93)

# read in community data file
data_file <- read_csv(file.path("sdl_veg", "data", "saddptqd.hh.data.csv"))

write.csv(data_file, file.path("sdl_veg", "figures", "mypretendfigures.csv"))

head(data_file)
dim(data_file) # 239864, 10

# look at data columns
colnames(data_file)
unique(data_file$LTER_site) # only NWT
unique(data_file$local_site) # only sad_grid
unique(data_file$year) #1989 - 2019 (17 yrs, but should be 16; 1996 should be 1995 for plot 37)
unique(data_file$plot) # 1-80 in a row, then also 101, 201, ... 801 for the A row (88 in total)
length(unique(data_file$point)) # I am not sure what this column represents - not the same as the length of the entire dataset - 8801
unique(data_file$x) # 5, 15,25, 35,45,55,65,75,85,95, NaN
unique(data_file$y) # same as x
unique(data_file$hit_type) # bottom, top, middle1, middle2, extra, middle3
unique(data_file$USDA_code) # 128
unique(data_file$USDA_name) # 132
# 2RF = rock fragments, 2LICHN = lichen, 2X = registration marker, 2LTR = Litter,  2BARE = bare ground, 2FORB = unknown forb, 2GRAM = unknown graminoid, 2MOSS = moss etc ... will need ti remove the unknown species as well as seperate out lichen, moss, rock, litter (everything starting with 2 is a non-species)
non_plant <- c("2RF", "2LICHN", "2X","2LTR","2BARE", "2HOLE", "2MOSS", "2SCATE")
unknown_species <- c("2FORB","2MOSS","2GRAM", "2UNKSC", "POA","2UNK", "CAREX", "2COMP") 

# look at patterns of hits across time
kinds_of_hits <- data_file %>% 
  group_by(plot, year) %>% 
  summarise(kinds_of_hits = unique(hit_type))
kinds_of_hits

# top hits
kinds_of_hits %>%
  filter(kinds_of_hits == "top") %>% 
  ggplot(aes(x = year, y = plot))+ 
  geom_point()+
  theme_classic()+
  ylim(0,80) # can remove this line to look at all plots, but just for ease
# there are some plots that never have top hits, there are also a small subset of plots that do not have tophits for only the first survey in 1989

#bottom hits
kinds_of_hits %>%
  filter(kinds_of_hits == "bottom") %>% 
  ggplot(aes(x = year, y = plot))+ 
  geom_point()+
  theme_classic()+
  ylim(0,80)
# all plots have bottom hits every year in the timeseries

# what are the plots that never have tophits? Are they just rocks? Yes
data_file %>% 
  group_by(plot) %>% 
  filter(all(hit_type == "bottom")) %>% 
  summarise(unique_species = unique(USDA_code)) %>% 
  as.data.frame()

# But what is going on in 1989? 
data_file %>% 
  filter(year == 1989) %>% 
  group_by(plot) %>% 
  filter(all(hit_type == "bottom")) %>% 
  summarise(unique_species = unique(USDA_code)) %>% 
  as.data.frame()

# are there more 'unknowns' in 1989? No - > 400 in 1989 and 1990, has decreased through time
data_file %>% 
filter(USDA_code %in% unknown_species) %>% 
group_by(year) %>%
tally() %>% 
  ggplot(aes(year, n))+
  geom_point()

# keep only bottom and top hits to allow for consistency across time
data_file <-
  data_file %>% 
  filter(hit_type == "bottom" | hit_type == "top")

# Additional subsetting --> top-most hit at each hit in each year (e.g. bottom hit if there is only a bottom hit, otherwise, top hit)
# checking...
summarized_hittype <-data_file %>% 
  arrange(year, plot, x, y, desc(hit_type)) %>% 
  group_by(year, plot, x,y) %>% 
  tally() %>% 
  arrange(n)

data_file %>% 
  filter(year == 1989, plot == 2, x == 5, y == 15)

data_file %>% 
  filter(year == 1989, plot == 1, x == 5, y == 5)

data_file <-data_file %>% 
  arrange(year, plot, x, y, desc(hit_type)) %>% 
  group_by(year, plot, x,y) %>% 
  slice(1)

# aggregate data_file to calculate the number of hits (adds top and bottom together)
number_hits <- data_file %>% 
  group_by(year, plot, USDA_code, USDA_name) %>% 
  summarise(hits = n())

# change year = 1996 to 1995 for plot 37 (typo in dataset)
number_hits[number_hits$year == 1996,"year"] <- 1995

# Fill in zeros for species/non-species not hit in a certain plot/year
data_with_zeros <- expand_grid(year = unique(number_hits$year), plot = unique(number_hits$plot), USDA_code = unique(number_hits$USDA_code))
data_with_zeros

abundance <- left_join(data_with_zeros, number_hits) %>% 
  dplyr::select(-USDA_name) 
abundance[is.na(abundance$hits), "hits"] <- 0
abundance

#### how has bare rock changed? ####
abundance %>% 
filter(USDA_code == "2RF") %>% 
ggplot(aes(x= year, y = hits))+
  geom_line(aes(group = plot), alpha = .1)+
  geom_smooth(method="glm.nb", se =F, color = "black")+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  ggtitle("Rock")

#abundance$year_factor <-as.factor(abundance$year) 
rock_temporal <- lmer(hits ~ year + (1|year) , data = abundance[abundance$USDA_code == "2RF",])
# won't converge

rock_temporal <- glmer.nb(hits ~ year + (1|year), data = abundance[abundance$USDA_code == "2RF",], control=glmerControl(optimizer="bobyqa"))
rock_temporal@optinfo[c("optimizer","control")]

summary(rock_temporal)# no significant change

#### how has total veg abundance changed? ####
total_abundance<-
abundance %>% 
  filter(!(USDA_code %in% non_plant)) %>% 
  group_by(year, plot) %>% 
  summarise(hits = sum(hits))

ggplot(total_abundance, aes(x= year, y = hits))+
  geom_point(alpha = .1, cex = 1)+
  geom_line(aes(group = plot), alpha = .1)+
  geom_smooth(method = "lm", color = "black", se = F)+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  #annotate(geom="text", size = 5, x=2010, y=100, label= "Estimate = 1.12, P-val < 00.1")+
  theme(text = element_text(size=15))+
  ggtitle("Total Vegetation Abundance")

veg_temporal <- lmer(hits ~ year + (1|year) + (1|plot), data = total_abundance)
fm2 <- lmer(hits ~ (1 | year) + (1|plot), data = total_abundance)
summary(veg_temporal) # total veg is significantly increasing!
AIC(veg_temporal)
anova(veg_temporal, fm2)

veg_resids <- as.data.frame(residuals(veg_temporal))
veg_resids$year <- total_abundance$year
plot(veg_resids$year, veg_resids$"residuals(veg_temporal)")
acf(veg_resids$"residuals(veg_temporal)")

#### how has bare ground changed? ####
abundance %>% 
filter(USDA_code == "2BARE") %>% 
ggplot(aes(x= year, y = hits))+
  geom_line(aes(group = plot), alpha = .1)+
  geom_smooth(method="glm.nb", se =F, color = "black")+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  #annotate(geom="text", size =5, x=2010, y=65, label= "Exp(Estimate) = 0.98, P-val < 0.001")+
  theme(text = element_text(size=15))+
  ggtitle("Bare Ground")

bare_temporal <- glmer.nb(hits ~ year + (1|plot), data = abundance[abundance$USDA_code == "2BARE",], control=glmerControl(optimizer="bobyqa"))
summary(bare_temporal) # bare ground is significantly decreasing!

#### how has lichen/mosses changed? ####
abundance %>% 
filter(USDA_code == "2LICHN") %>% 
ggplot(aes(x= year, y = hits))+
  geom_line(aes(group = plot), alpha = .1)+
  geom_smooth(method = "glm.nb", color = "black", se = F)+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  #annotate(geom="text", x=2010, y=80, size =5, label= "Exp(Estimate) = 0.99, P-val = 0.01")+
  theme(text = element_text(size=15))+
  ggtitle("Lichen")

library(GLMMadaptive)
lichen_temporal <- mixed_model(fixed = hits ~ year, random = ~ 1 | year, data = abundance[abundance$USDA_code == "2LICHN",], family = negative.binomial())
summary(lichen_temporal) # lichen is significantly decreasing!

abundance %>% 
filter(USDA_code == "2MOSS") %>% 
ggplot(aes(x= year, y = hits))+
  geom_line(aes(group = plot), alpha = .1)+
  #geom_smooth(method = "glm.nb",se=F, color = "black")+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  theme(text = element_text(size=15))+
  ggtitle("Moss")

moss_temporal <- glm.nb(hits ~ year, data = abundance[abundance$USDA_code == "2MOSS",])
summary(moss_temporal)

### how have shrubs changed ####
# shrubs in the saddle 
shrubs <- c("SANI8", "SAPL2", "SAPE18", "SAGL")

abundance %>% 
filter(USDA_code %in% shrubs) %>% 
  group_by(year, plot) %>% 
  summarise(hits = sum(hits)) %>% 
  ggplot(aes(x= year, y = hits))+
  geom_line(aes(group = plot), alpha = .5)+
  #geom_smooth(method = "glm.nb", color = "black", se = F)+
  xlab("Year")+
  ylab("Abundance")+
  theme_classic()+
  #annotate(geom="text", x=2010, y=80, size =5, label= "Exp(Estimate) = 0.99, P-val = 0.01")+
  theme(text = element_text(size=15))+
  ggtitle("Shrubs")

shrubs_abundance <- 
  abundance %>% 
 filter(USDA_code %in% shrubs) %>% 
  group_by(year, plot) %>% 
  summarise(hits = sum(hits)) 

shrub_temporal <- glm.nb(hits ~ year, data = shrubs_abundance)
summary(shrub_temporal)

#### changes in richness ####
richness <-
  abundance %>% 
  filter(!(USDA_code %in% non_plant)) %>% 
  filter(!(USDA_code %in% unknown_species)) %>%
  filter(hits > 0) %>% 
  group_by(year, plot) %>% 
  summarise(richness = length(unique(USDA_code)))
richness

ggplot(richness, aes(x= year, y = richness))+
  geom_line(aes(group = plot), alpha = .1)+
  geom_smooth(method = "lm",se=F, color = "black")+
  xlab("Year")+
  ylab("Richness")+
  theme_classic()+
  #annotate(geom="text", x=2010, y=25,  size =5, label= "Estimate = 0.09, P-val < 0.001")+
  theme(text = element_text(size=15))+
  ggtitle("Richness")

richness_model <- lmer(richness ~ year + (1|plot), data = richness)
summary(richness_model)
fm2 <- lmer(richness ~ (1 | year), data = richness)
anova(richness_model, fm2) 

#### Communities through time ####
# need to make year/plot specific rows
community<- 
  abundance %>% 
  filter(!(USDA_code %in% non_plant)) %>% 
  filter(!(USDA_code %in% unknown_species)) %>% 
  filter(!(USDA_code %in% "CASCS2")) # there are duplicates of this species
community

subset_community<-
community %>% 
  group_by(year, plot) %>% 
  mutate(totals = sum(hits)) %>% 
  filter(totals > 0)

# remove 1989 due to methodological issues? Nope
#subset_community <-
#  subset_community %>% 
#  filter(year >1989)

# multivariate change
?multivariate_change
mumti_change <- multivariate_change(df=subset_community, time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
mumti_change 

plot1<- mumti_change %>% 
ggplot(aes(year2, composition_change)) +
geom_point()+
geom_line()+
theme_classic()+
xlab("Year")+
ylab("Compositional Change Relative to 1989")+
ggtitle("Saddle Veg Temporal Change")+
ylim(0, .14)  

plot2<- mumti_change %>% 
ggplot(aes(year2, dispersion_change)) +
geom_point()+
geom_line()+
theme_classic()+
xlab("Year")+
ylab("Dispersion Change Relative to 1989")+
#ylim(-.025, 0.01)+
geom_hline(yintercept = 0)+
ggtitle("Saddle Veg Spatial Change")

plot1
plot2

cowplot::plot_grid(plot1, plot2,nrow =2)

summary(lm(composition_change~year2,data=mumti_change))
summary(lm(dispersion_change~year2,data=mumti_change))


# merge in info about veg classifications
# FF=fellfield, DM=dry meadow, MM=moist meadow, ST=shrub tundra, SB=snowbed, WM=wet meadow, SF=snowfence
#veg_classes <- read_csv("Saddle_Analyses/data/raw/saddgrid_npp.hh.data.csv")
veg_classes <- read_csv(here("data", "saddgrid_npp.hh.data.csv")) 
unique(veg_classes$year)
unique(veg_classes$grid_pt)
unique(veg_classes$veg_class)
#why do some plots have an "A" in them? 

#subset to 
veg_classes <-
  veg_classes %>% 
  filter(year == 2010, subsample == "A")
veg_classes<-veg_classes[,c(5,7)]
colnames(veg_classes) <- c("plot", "veg_class")

# merge community data and veg classes
subset_community$plot<-as.character(subset_community$plot) #needed because could not join with plot as a double
veg_classes$plot<- as.character(veg_classes$plot)
community_class <- left_join(subset_community,veg_classes)
community_class
unique(community_class$veg_class)

#MM
MM_change<- 
community_class %>% 
  filter(veg_class == "MM") %>% 
multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
MM_change$class <- "Moist Meadow"
MM_change

#SB
SB_change<-
community_class %>% 
  filter(veg_class == "SB") %>% 
multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
SB_change$class <- "Snow Bed"
SB_change

#WM
WM_change<- 
community_class %>% 
  filter(veg_class == "WM") %>%   
multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
WM_change$class <- "Wet Meadow"
WM_change

#SF
SF_change<- 
community_class %>% 
  filter(veg_class == "SF") %>% 
multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
SF_change$class <- "Snow Fence"
SF_change

#DM
DM_change<-
community_class %>%
 filter(veg_class == "DM") %>%
 multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
DM_change$class <- "Dry Meadow"
DM_change

#FF 
FF_change<-
 community_class %>%
 filter(veg_class == "FF") %>%
 multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
FF_change$class <- "Fellfield" 
FF_change

#ST
ST_change<-
community_class %>%
 filter(veg_class == "ST") %>%
 multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = 1989)
ST_change$class <-"Shrub Tundra"
ST_change

all_veg_classes <- rbind(MM_change, SB_change, WM_change, DM_change, FF_change)
all_veg_classes

plot3<-all_veg_classes %>% 
ggplot(aes(year2, composition_change, color = class)) +
geom_point()+
geom_line(lwd = 2)+
theme_cowplot()+
xlab("Year")+
ylim(0,.38)+
ylab("Compositional Change Relative to 1989")+
      scale_colour_viridis_d(option = "turbo", labels = c("Dry Meadow", "FellField", "Moist Meadow", "Snowbed", "Wet Meadow"))+
  theme(legend.title=element_blank(), text = element_text(size=15))+
ggtitle("Saddle Vegetation Temporal Change")
plot3

plot4<-all_veg_classes %>% 
ggplot(aes(year2, dispersion_change, color = class)) +
geom_point()+
geom_line(lwd = 1.5)+
theme_cowplot()+
xlab("Year")+
ylab("Dispersion Relative to 1989")+
geom_hline(yintercept = 0)+
ggtitle("Saddle Vegetation Spatial Change")+
  scale_colour_viridis_d(option = "turbo", labels = c("Dry Meadow", "FellField", "Moist Meadow", "Snowbed", "Wet Meadow"))+
  theme(legend.title=element_blank(), text = element_text(size=20))
plot4
#cowplot::plot_grid(plot3, plot4,nrow =2)

# Look saddle temporal change across climate variables
all_veg_classes

# bring in snow melt dates
snow_melt_dates <- read_csv("Saddle_Analyses/data_deriv/snowmelt_est.csv")
snow_melt_dates <- snow_melt_dates[,c(1,2,6)]
colnames(snow_melt_dates)[2] <- "plot"
colnames(snow_melt_dates)[3] <- "melt_date"
snow_melt_dates

# calculate a veg-class average melt date through time
snow_melt_dates <- left_join(snow_melt_dates, veg_classes)
snow_melt_dates$class <- "Place_Holder"
snow_melt_dates[which(snow_melt_dates$veg_class == "MM"), "class" ] <- "Moist Meadow"
snow_melt_dates[which(snow_melt_dates$veg_class == "FF"), "class" ] <- "Fellfield"
snow_melt_dates[which(snow_melt_dates$veg_class == "DM"), "class" ] <- "Dry Meadow"
snow_melt_dates[which(snow_melt_dates$veg_class == "WM"), "class" ] <- "Wet Meadow"
snow_melt_dates[which(snow_melt_dates$veg_class == "SB"), "class" ] <- "Snow Bed"
snow_melt_dates

# make an annual melt date dependent on veg type or not
snow_melt_veg <- 
  snow_melt_dates %>% 
  group_by(veg_class, class, year) %>% 
  summarise(melt_veg_specific = mean(melt_date))
snow_melt_veg 
colnames(snow_melt_veg)[3] <- "year2"

snow_melt_annual <- 
  snow_melt_dates %>% 
  group_by(year) %>% 
  summarise(melt_annual = mean(melt_date))
snow_melt_annual 
colnames(snow_melt_annual)[1] <- "year2"


veg_snow <- left_join(all_veg_classes, snow_melt_veg)
veg_snow <- left_join(veg_snow, snow_melt_annual)

veg_snow %>% 
ggplot(aes(melt_annual, composition_change, color = class)) +
geom_point()+
geom_smooth(data = veg_snow[veg_snow$class == "Wet Meadow",], method = "lm", se = F)+
geom_smooth(data = veg_snow[veg_snow$class != "Wet Meadow",], method = "lm", se = F, lty = "dashed")+
theme_cowplot()+
xlab("Annual Average Melt Date")+
ylab("Compositional Change Relative to 1989")

fit_melt <- lmer(composition_change ~ melt_annual*class + (1|year2), data = veg_snow)
summary(fit_melt)
emtrends(cnm_GDD_veg_model, specs = "veg_class", var = "GDD")

veg_snow %>% 
ggplot(aes(melt_veg_specific, composition_change, color = class)) +
geom_point()+
geom_smooth(data = veg_snow[veg_snow$class == "Wet Meadow",], method = "lm", se = F)+
geom_smooth(data = veg_snow[veg_snow$class != "Wet Meadow",], method = "lm", se = F, lty = "dashed")+
theme_cowplot()+
xlab(" Veg-class Specific Annual Average Melt Date")+
ylab("Compositional Change Relative to 1989")

summary(lmer(composition_change ~ melt_veg_specific*class + (1|year2), data = veg_snow))

#bring in temperature data
temp <- read_csv("c1_d1_sdl_clim/homogenize_climdat/data/sdl_temp_1981-2020_draft.csv")
#temp$date <-ymd(temp$date); temp$year <- year(temp$date); temp$month <- month(temp$date)

# calculate GDD
GDD <- 
  temp %>% 
  filter(yr >= 1989, metric == 'mean') %>% 
  filter(adjusted_airtemp > 0) %>% 
  group_by(yr) %>% 
  summarise(GDD = sum(adjusted_airtemp))
colnames(GDD)[1] <- "year2"
veg_temp <- left_join(all_veg_classes, GDD)
veg_temp

to_use_col <- turbo(5, alpha = 1, begin = 0, end = 1, direction = 1)
 

VegChange_GDD<- 
veg_temp %>% 
ggplot(aes(GDD, composition_change, color = class)) +
geom_point(size = 2)+
geom_smooth(data = veg_temp[!veg_temp$class %in% c("Dry Meadow","Moist Meadow"),], method = "lm", se = F)+
geom_smooth(data = veg_temp[veg_temp$class %in% c("Dry Meadow","Moist Meadow"),], method = "lm", se = F, lty = "dashed")+
theme_classic()+
xlab("GDD")+
ylab("Compositional change")+
  scale_color_viridis_d(option = "turbo", labels = c("DM", "FF", "MM", "SB", "WM"))+
 #scale_fill_viridis_d(option = "turbo")+
theme(legend.position = "bottom", legend.text.align = 0, text = element_text(size=16), legend.title = element_blank(),)+
  annotate("text", x = 1360, y = .12, label = "DM", color = to_use_col[1], size = 3)+
  annotate("text", x = 1360, y = .16, label = "MM", color = to_use_col[3], size = 3)+
  annotate("text", x = 1360, y = .21, label = "SB", color = to_use_col[4], size = 3)+
  annotate("text", x = 1360, y = .3, label = "FF", color = to_use_col[2], size = 3)+
  annotate("text", x = 1360, y = .33, label = "WM", color = to_use_col[5], size = 3)+
  theme(legend.position = "none", legend.text.align = 0, axis.text.x= element_text(color = "black"), axis.text.y= element_text(color = "black"))
#+ggtitle("B")

VegChange_GDD


ggsave(plot = VegChange_GDD ,filename = "Saddle_Analyses/figures/CompChange_GDD.jpg", width = 3, height= 3, units = "in")
  

fit_GDD <- lmer(composition_change ~ GDD*class + (1|year2), data = veg_temp)
summary(fit_GDD)
emtrends(fit_GDD , specs = "class", var = "GDD")

# consider melt date corrected GDD
GDD_melt_corrected <-
temp %>% 
  filter (yr >= 1992) %>%  #filter to later years when snow melt available
  mutate(date_ord = yday(date)) %>% 
  left_join(., snow_melt_veg, by = c("yr" = "year2")) %>% 
  filter(class != "Place_Holder") %>% 
  filter(date_ord > melt_veg_specific) %>%  # remove all days prior to the snow melt date 
  filter(adjusted_airtemp > 0) %>%  # remove all days below zero
  group_by(yr, class) %>% 
  summarise(GDD_snow = sum(adjusted_airtemp))
GDD_melt_corrected
colnames(GDD_melt_corrected)[1] <- "year2"

veg_temp_corrected<- left_join(all_veg_classes, GDD_melt_corrected)
veg_temp_corrected

veg_temp_corrected %>% 
ggplot(aes(GDD_snow, composition_change, color = class)) +
geom_point()+
geom_smooth(data = veg_temp_corrected[!veg_temp_corrected$class %in% c("Dry Meadow", "Moist Meadow", "Snow Bed"),], method = "lm", se = F)+
geom_smooth(data = veg_temp_corrected[veg_temp_corrected$class %in% c("Dry Meadow", "Moist Meadow", "Snow Bed"),], method = "lm", se = F, lty = "dashed")+
theme_cowplot()+
xlab("Growing Degree-Days (melt_corrected)")+
ylab("Compositional Change Relative to 1989")+
    scale_color_viridis_d(option = "turbo")

fit_melt_corrected <- lmer(composition_change ~ GDD_snow*class + (1|year2), data = veg_temp_corrected)
summary(fit_melt_corrected)
emtrends(fit_melt_corrected , specs = "class", var = "GDD_snow")

# Examine how GDD and GGD_corrected relate to changes in total vegetation abundance, bare ground, and lichens/moss
#total abundance
total_abundance<-
abundance %>% 
  filter(!(USDA_code %in% non_plant)) %>% 
  group_by(year, plot) %>% 
  summarise(hits = sum(hits))
colnames(total_abundance)[1] <- "year2"

# make a GDD_melt not specific to veg classes
GDD_melt_corrected <-
temp %>% 
  filter (year > 1992) %>%  #filter to later years when snow melt available
  mutate(date_ord = yday(date)) %>% 
  left_join(., snow_melt_annual , by = c("year" = "year2")) %>% 
  filter(date_ord > melt_annual) %>%  # remove all days prior to the snow melt date 
  filter(average > 0) %>%  # remove all days below zero
  group_by(year) %>% 
  summarise(GDD_snow = sum(average))
GDD_melt_corrected
colnames(GDD_melt_corrected)[1] <- "year2"

total_veg_clim <-
total_abundance %>% 
  left_join(., GDD) %>% 
  left_join(., GDD_melt_corrected) %>% 
  left_join(., snow_melt_annual) 
total_veg_clim

total_veg_clim %>% 
filter(hits > 0) %>% 
ggplot(aes(GDD, hits)) +
geom_point(alpha = .1)+
geom_smooth(method = "lm", color = "black", se = F)+
theme_cowplot()+
xlab("Growing Degree-Days")+
ylab("Total Vegetation Abundance")

veg_GDD <- lmer(hits ~ GDD + (1|year2) + (1|plot), data = total_veg_clim)
summary(veg_GDD)

total_veg_clim %>% 
filter(hits > 0) %>% 
ggplot(aes(GDD_snow, hits)) +
geom_point(alpha = .1)+
geom_smooth(method = "lm", color = "black", se = F)+
theme_cowplot()+
xlab("Growing Degree-Days (melt_corrected)")+
ylab("Total Vegetation Abundance")

veg_GDD_melt <- lmer(hits ~ GDD_snow + (1|year2) + (1|plot), data = total_veg_clim)
summary(veg_GDD_melt)

total_veg_clim %>% 
filter(hits > 0) %>% 
ggplot(aes(melt_annual, hits)) +
geom_point(alpha = .1)+
geom_smooth(method = "lm", color = "black", se = F,lty = "dashed")+
theme_cowplot()+
xlab("Annual Melt Date")+
ylab("Total Vegetation Abundance")

veg_melt <- lmer(hits ~ melt_annual + (1|year2) + (1|plot), data = total_veg_clim)
summary(veg_melt)

#####  look at lichens&mosses ####
lichen_mosses<-
abundance %>% 
  filter(USDA_code %in% c("2LICHN", "2MOSS")) %>% 
  pivot_wider(names_from = USDA_code, values_from = hits)
colnames(lichen_mosses)<- c("year2", "plot", "lichen", "moss")
lichen_mosses

lichen_clim <- left_join(lichen_mosses, total_veg_clim) 
lichen_clim

lichen_clim %>% 
ggplot() +
geom_point(aes(GDD, lichen), alpha = .5, color = "lightgreen")+
geom_point(aes(GDD, moss), alpha = .5, color = "forestgreen")+
theme_cowplot()+
geom_smooth(aes(GDD, lichen), color = "lightgreen", method = "lm", color = "black", se = F)+
geom_smooth(aes(GDD, moss), color = "forestgreen", method = "lm", color = "black", se = F)+
xlab("Growing Degree-Days")+
ylab("Moss & Lichen Abundance")

lichen_GDD <- lmer(lichen ~ GDD + (1|year2) + (1|plot), data = lichen_clim)
summary(lichen_GDD)
moss_GDD <- lmer(moss ~ GDD + (1|year2) + (1|plot), data = lichen_clim)
summary(moss_GDD)

lichen_clim %>% 
ggplot() +
geom_point(aes(GDD_snow, lichen), alpha = .5, color = "lightgreen")+
geom_point(aes(GDD_snow, moss), alpha = .5, color = "forestgreen")+
theme_cowplot()+
geom_smooth(aes(GDD_snow, lichen), color = "lightgreen", method = "lm", color = "black", se = F)+
geom_smooth(aes(GDD_snow, moss), color = "forestgreen", method = "lm", color = "black", se = F)+
xlab("Growing Degree-Days (melt_corrected)")+
ylab("Moss & Lichen Abundance")

lichen_GDD_melt <- lmer(lichen ~ GDD_snow + (1|year2) + (1|plot), data = lichen_clim)
summary(lichen_GDD_melt)
moss_GDD_melt<- lmer(moss ~ GDD_snow + (1|year2) + (1|plot), data = lichen_clim)
summary(moss_GDD_melt)

lichen_clim %>% 
ggplot() +
geom_point(aes(melt_annual, lichen), alpha = .5, color = "lightgreen")+
geom_point(aes(melt_annual, moss), alpha = .5, color = "forestgreen")+
theme_cowplot()+
geom_smooth(aes(melt_annual, lichen), color = "lightgreen", method = "lm", color = "black", se = F, lty = "dashed")+
geom_smooth(aes(melt_annual, moss), color = "forestgreen", method = "lm", color = "black", se = F, lty = "dashed")+
xlab("Annual Melt Date")+
ylab("Moss & Lichen Abundance")

lichen_melt <- lmer(lichen ~ melt_annual + (1|year2) + (1|plot), data = lichen_clim)
summary(lichen_melt)
moss_melt<- lmer(moss ~ melt_annual + (1|year2) + (1|plot), data = lichen_clim)
summary(moss_melt)


#### Loop for different veg types relative to all the years #####
community_class
classes <- unique(community_class$veg_class); classes
years <- unique(community_class$year); years

# how many plots per veg type?
community_class %>% 
  filter(year == 2010) %>% 
  group_by(veg_class) %>% 
  summarise(num_plots = length(unique(plot)))
# remove shrub tundra and snow fence as only 1  and 2 plots respectively
community_classes_subset <- community_class %>% 
  filter(veg_class != "SF") %>% 
  filter(veg_class != "ST")

# make a loop of different reference years
relative_to_what <-data.frame()

for(i in 1:length(classes)){
  for(j in 1:length(years)){
   
section<- 
  community_classes_subset %>% 
  filter(veg_class == classes[i]) %>% 
  multivariate_change(time.var = "year", species.var = "USDA_code",abundance.var = "hits", replicate.var = "plot", reference.time = years[j])
  section$class <- classes[i]
  
relative_to_what <- rbind(relative_to_what, section)
    
  }
}
  
head(relative_to_what)
relative_to_what$diff <- abs(relative_to_what$year - relative_to_what$year2)

relative_to_what %>% 
  ggplot(aes(diff, composition_change)) + 
        geom_point()+
      geom_smooth(se = F) + 
  facet_grid(.~class)+
  xlab("Difference in Years")+
  theme_bw()

relative_to_what %>% 
  ggplot(aes(diff, dispersion_change)) + 
        geom_point()+
  facet_grid(.~class)+
  theme_bw()


###############################
#EXTRACTING ORDINATION COORDS FOR ALL PLOT, YEAR COMBINATIONS
multivariate_change
pca_centers
split_apply_combine
# hmmm function from codyn may not be of much help
classes <- unique(community_class$veg_class); classes
all_centroids <- data.frame()

for(i in 1:length(classes)){
# subset by veg class
test <- community_class %>% 
  filter(veg_class == classes[i]) %>% 
  unite(year_plot, year, plot, sep = "_")
test

test_matrix <- sample2matrix(test[,c(1,3,2)])
test_dis <- vegdist(test_matrix, method = "bray")

# ordinate with PCoA (Bray-Curtis dissimilarity) data for all years, plots in that veg class
test_PCoA <- pcoa(test_dis)
barplot(test_PCoA$values$Relative_eig[1:10])
biplot.pcoa(test_PCoA)

# record locations of all plots for each year
plot_coords <- tibble(id = rownames(test_matrix), A1 =  test_PCoA$vectors[,1], A2 =  test_PCoA$vectors[,2])

plot_coords <-
  plot_coords %>% 
  separate(id, into = c("year", "plot"))
  
# find centroid for each year (group mean) and record 
centroids <- 
  plot_coords %>% 
  group_by(year) %>% 
  summarise(C1 = mean(A1), C2 = mean(A2))

centroids$veg_class <- classes[i]

all_centroids<- rbind(all_centroids, centroids)

}


all_centroids %>% 
  filter(veg_class != "SF") %>% 
  ggplot(aes(C1, C2))+
  geom_text(aes(label = year))+
  xlab("PC1")+
  ylab("PC2")+
  theme_classic()+
  geom_link2(lineend = 'round', n = 500, color = "gray")+
  facet_wrap(vars(veg_class), scales = "free")


# make ordi plot for all veg classes together
all <- community_class %>% 
  filter(veg_class != "SF") %>% 
  unite(year_plot, year, plot, sep = "_")
all

all_matrix <- sample2matrix(all[,c(1,3,2)])
all_dis <- vegdist(all_matrix, method = "bray")

# ordinate with PCoA (Bray-Curtis dissimilarity) data for all years, plots in that veg class
all_PCoA <- pcoa(all_dis)
barplot(all_PCoA$values$Relative_eig[1:10])
biplot.pcoa(all_PCoA)

# record locations of all plots for each year
plot_coords_all <- tibble(id = rownames(all_matrix), A1 =  all_PCoA$vectors[,1], A2 =  all_PCoA$vectors[,2])

plot_coords_all <-
  plot_coords_all %>% 
  separate(id, into = c("year", "plot"))
  
# find centroid for each year (group mean) and record 
centroids_allVeg <- 
  plot_coords_all %>% 
  group_by(year) %>% 
  summarise(C1 = mean(A1), C2 = mean(A2))

centroids_allVeg %>% 
  ggplot(aes(C1, C2))+
  geom_text(aes(label = year))+
  xlab("PC1")+
  ylab("PC2")+
  theme_classic()+
  geom_link2(lineend = 'round', n = 500, color = "gray")

##################### Old code ####

# #overall community turnover
# saddle_turnover <- turnover(df = subset_community, time.var = "year", species.var = "USDA_code", abundance.var = "hits", replicate.var = "plot")
# saddle_turnover
# 
# saddle_turnover %>% 
#   group_by(year) %>% 
#   summarise(sum.total = sum(total)) %>% 
# ggplot(aes(x = year, y = sum.total))+
#   geom_line(lwd = 1.5)+
#   theme_classic()+
#     xlab("Year")+
#   ylab("Average Turnover")+
#   theme(text = element_text(size=15))+
#   ggtitle("Community Turnover")
# 
# saddle_turnover %>% 
# ggplot(aes(x = year, y = total, group = plot))+
#   geom_line(lwd=.2)+
#   theme_classic()+
#     xlab("Year")+
#   ylab("Average Turnover")+
#   theme(text = element_text(size=15))+
#   ggtitle("Community Turnover")
# 
# 
# # disappearances
# ?turnover
# saddle_dis <- turnover(df = subset_community, time.var = "year", species.var = "USDA_code", abundance.var = "hits", replicate.var = "plot", metric = "disappearance")
# saddle_dis
# 
# saddle_dis %>% 
#   #group_by(year) %>% 
#   #summarise(sum.total = sum(disappearance)) %>% 
# ggplot(aes(x = year, y = disappearance, group = plot))+
#   geom_line()+
#   theme_classic()+
#     xlab("Year")+
#   ylab("Average Turnover")+
#   theme(text = element_text(size=15))+
#   #ylim(5,18)+
#   ggtitle("Disappearances")
# 
# # appearances
# saddle_app <- turnover(df = subset_community, time.var = "year", species.var = "USDA_code", abundance.var = "hits", replicate.var = "plot", metric = "appearance")
# saddle_app
# 
# saddle_app %>% 
#   group_by(year) %>% 
#   summarise(sum.total = sum(appearance)) %>% 
# ggplot(aes(x = year, y = sum.total))+
#   geom_line(lwd = 1.5)+
#   theme_classic()+
#     xlab("Year")+
#   ylab("Average Turnover")+
#   theme(text = element_text(size=15))+
#   ylim(5,18)+
#   ggtitle("Appearances")
# 
# # rank shifts
# saddle_shifts <- rank_shift(df = subset_community, time.var = "year", species.var = "USDA_code", abundance.var = "hits", replicate.var = "plot")
# saddle_shifts$year <- as.numeric(substr(saddle_shifts$year_pair, 6, 9))
# 
# saddle_shifts %>% 
#   group_by(year) %>% 
#   summarise(sum.MRS = sum(MRS)) %>% 
# ggplot(aes(x = year, y = sum.MRS))+
#   geom_line(lwd = 1.5)+
#   theme_classic()+
#     xlab("Year")+
#   ylab("MRS")+
#   theme(text = element_text(size=15))+
#   ggtitle("Rank Shifts")
# 
# 
# saddle_shifts %>% 
# ggplot(aes(x = year, y = MRS, group = plot))+
#   geom_line(lwd = .2)+
#   theme_classic()+
#     xlab("Year")+
#   ylab("MRS")+
#   theme(text = element_text(size=15))+
#   ggtitle("Rank Shifts")
# 
# 
# 
# 
# 
# 
