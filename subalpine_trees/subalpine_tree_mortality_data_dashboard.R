# ============================================================================
# Niwot Ridge LTER Data Dashboard: Subalpine Forest Tree Mortality Figures
# ============================================================================
#
# Purpose:
#   Reads permanent forest plot census data (tree-level mortality records
#   across repeated surveys, 1982-2019), calculates annual percent
#   mortality by species and topographic moisture class, 
#   and produces two figures:
#
#     1. subalpine_tree_mortality.jpg          
#                 % annual mortality by census period and topographic moisture
#                 class, for Subalpine Fir (ABLA)and Engelmann Spruce (PIEN)
#     2. subalpine_tree_mortality_long-term_mean.jpg   
#                 Same comparison, faceted by species and topographic class,
#                 colored by direction relative to each species/class's 
#                 long-term mean mortality rate
#
# Input:
#   data/*.csv  -- permanent forest plot data (EDI knb-lter-nwt.207): plot
#                  metadata, tree-level census records, and seedling counts.

#
# Output:
#   Two JPG figures written to plot_dir (see file names above).
#
# Notes:
#   Restricted, per R. Andrus's recommendation, to plots BW3, MRS4, MRS5,
#   BL6, and MRS7, and to Subalpine Fir (ABLA) and Engelmann Spruce (PIEN)
#   -- Quaking Aspen (POTR) has too few trees for confidence, and Lodgepole
#   Pine (PICO) mortality appears largely density-dependent rather than
#   climate-driven.
#
# TODO
#   - double check the "one alive at end of period" vs "died during period"
#     mortality accounting logic against the source publication
#   - revisit whether the 100%-mortality plots (single last tree died) are
#     handled the way you want -- see mortrates_all / PercAnnMort_adj
#
# Original Code and Analysis: Sarah C. Elmendorf(2023), building on an
#   earlier tree-mortality-by-PC1 workflow by Cliff Bueno de Mesquita, and
#   on download/subsetting utilities by C. Bueno de Mesquita
# Revised and Updated for the Data Dashboard by Anne Marie Panetta (2026)

# Data source: Veblen, T., R. Andrus, and R. Chai. 2021. Permanent forest
#   plot data from 1982-2019 at Niwot Ridge, Colorado ver 5. Environmental
#   Data Initiative. 
#   https://doi.org/10.6073/pasta/48fa11a8f5bc6541b0472bc3fd4c0c71
# ============================================================================

# -- SETUP --------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(EDIutils)  # tools for interacting with EDI's data package API
library(plyr)       # ddply(), used inside summarySE()

options(stringsAsFactors = FALSE)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_dir <- "data"
plot_dir <- "plots"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Set to TRUE to (re)download source data from EDI. Only needs to be run
# once, or when a newer package version is released.
download_data <- FALSE



# -- HELPER FUNCTION ------------------------------------------------------
# Summarizes a numeric variable by group: N, mean, sd, standard error, and
# a confidence interval, in the style of the Cookbook for R's summarySE().

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, returns a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#TODO 


# -- DOWNLOAD DATA FROM EDI ---------------------------------------------
# Note: if you have already downloaded SOME data, read_data_package_archive()
# will not overwrite existing files. Clear the /data directory first if you
# need to force a fresh download.
if (download_data) {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  scope <- "knb-lter-nwt"  # Niwot scope
  
  # 207 -- permanent forest plot data (plot metadata, tree census, and
  # seedling tables), with current citation (update as package version
  # changes):
  #   Veblen, T., R. Andrus, and R. Chai. 2021. Permanent forest plot data
  #   from 1982-2019 at Niwot Ridge, Colorado ver 5. Environmental Data
  #   Initiative. https://doi.org/10.6073/pasta/48fa11a8f5bc6541b0472bc3fd4c0c71
  for (id in c("207")) {
    revision <- list_data_package_revisions(scope, id, filter = "newest")
    packageID <- paste(scope, id, revision, sep = ".")
    
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID)) # confirm the version used
  }
  
  for (fname in list.files(data_dir, pattern = "knb-lter.*zip", 
                           full.names = TRUE)) {
    unzip(zipfile = fname, exdir = data_dir)
  }
  
  # After the first download, list what actually landed in data_dir and
  # update the three file names in READ DATA below to match.
  print(list.files(data_dir, pattern = "\\.csv$", ignore.case = TRUE))
}
# -- READ DATA ------------------------------------------------------------
# TODO: the package includes three tables (plot metadata, tree-level
# census records, and seedling counts). Confirm the actual file names in
# data_dir (see the list.files() call above, or run it manually) and
# update the three paths below -- these are placeholders.

PP_plot_data <- read_csv(file.path(data_dir, "PP_plot_data.tv.data.csv"), 
                         na = na_vals)
PP_tree_data <- read_csv(file.path(data_dir, "PP_tree_data.tv.data.csv"), 
                         na = na_vals)
PP_seedling_data <- read_csv(file.path(data_dir, 
                            "PP_seedling_data.tv.data.csv"), 
                             na = na_vals)

# -- PREP PLOT METADATA -----------------------------------------------------

# Topographic position and stand age for each plot / plot group, from
# https://www.sciencedirect.com/science/article/pii/S0378112714007476#t0005
#
#   Gap plots 6-10 and 12-21 are over a broad area west of near MRS4
#     (old, xeric/mesic)
#   Gap plots 11 and 22-30 are near BL6 (old, mesic)
#   MRS8, 9 & 10 (xeric/mesic, young)
#   Gap plots 1-5 are in a tight cluster on the north side of the ski
#     trail north of MRS1 (young, xeric)
extra_plot_info <- data.frame(
  plot = c("BW2", "BW3", "MRS4", "MRS5", "BL6", "MRS7", "MRS4_gap", "BL6_gap", 
           "MRS1", "MRS8-9-10", "MRS1_gap"),
  topo = c("Xeric", "Mesic", "Xeric", "Hydric", "Mesic", "Xeric", "Xeric/Mesic", 
           "Mesic", "Xeric", "Xeric/Mesic", "Xeric"),
  age = c("old", "old", "old", "old", "old", "old", "old", "young", "young", 
          "young", "young")
)

# Individual "gap" sub-plots are pooled into their parent plot group.
PP_tree_data <- PP_tree_data %>%
  dplyr::mutate(plot_group = plot) %>%
  dplyr::mutate(plot_group = replace(plot_group, plot_group %in% c("Gap1", 
                            "Gap2", "Gap3", "Gap4", "Gap5"), "MRS1_gap")) %>%
  dplyr::mutate(plot_group = replace(plot_group, plot_group %in% 
                              c("MRS8", "MRS9", "MRS10"), "MRS8-9-10")) %>%
  dplyr::mutate(plot_group = replace(plot_group, plot_group %in% c(
    "Gap6", "Gap7", "Gap8", "Gap9", "Gap10",
    "Gap12", "Gap13", "Gap14", "Gap15", "Gap16", 
    "Gap17", "Gap18", "Gap19", "Gap20", "Gap21"
  ), "MRS4_gap")) %>%
  dplyr::mutate(plot_group = replace(plot_group, plot_group %in% c(
    "Gap11", "Gap22", "Gap23", "Gap24", "Gap25", "Gap26", "Gap27", 
    "Gap28", "Gap29", "Gap30"
  ), "BL6_gap")) %>%
  dplyr::left_join(., PP_plot_data %>% dplyr::select(plot, install_year))

# Fill in install years for plot groups not resolved above (from the
# source publication).
PP_tree_data$install_year[PP_tree_data$plot_group == "MRS1_gap"] <- 1983
PP_tree_data$install_year[PP_tree_data$plot_group == "BL6_gap"] <- 1983
PP_tree_data$install_year[PP_tree_data$plot_group == "MRS4_gap"] <- 1983

# For trees with no recorded start of their "dead" period (dp_start), but
# with a recorded end (dp_end) and marked dead, assume the dead period
# began at plot installation.
PP_tree_data$dp_start[is.na(PP_tree_data$dp_start) & 
                        PP_tree_data$dead == 1 & !is.na(PP_tree_data$dp_end)] <-
  PP_tree_data$install_year[is.na(PP_tree_data$dp_start) 
                        & PP_tree_data$dead == 1 & !is.na(PP_tree_data$dp_end)]

PP_tree_data$duration <- (PP_tree_data$dp_end - PP_tree_data$dp_start) + 1

PP_tree_data$dbhmean <- rowMeans(PP_tree_data[c("dbh1", "dbh3")], na.rm = TRUE)


# -- ASSIGN CENSUS-PERIOD STATUS -------------------------------------------
# For each tree at each census period, assign a status:
#   a  = alive during this period
#   d  = died during this period
#   dp = died in a previous period (still absent/dead in this one)
#
# Note on mortality accounting: the published data calculated percent
# mortality based on the number of trees at plot installation, but the
# more defensible denominator is the number of trees alive at the start
# of each period -- i.e. mortality = dead / (dead + alive), not
# dead / alive. This is corrected further down.

PP_tree_data <- PP_tree_data %>% tidyr::unite(., period, dp_start, dp_end, 
                                              sep = "-")

names(PP_tree_data)[2] <- c("tree.")  # rename the tree ID column

mortality <- tidyr::expand(PP_tree_data, tidyr::nesting(plot, tree.), period) %>%
  dplyr::inner_join(., PP_tree_data %>% select(plot, period) %>% distinct())  
# drop census periods that don't actually exist for a given plot

mortality <- mortality %>%
  dplyr::left_join(., PP_tree_data %>% dplyr::select(
    plot, tree., spec, period, dead, hc1, biotic_attack, biotic_agent, plot_group
  ) %>% dplyr::rename(mort_period = period)) %>%
  dplyr::arrange(plot, tree., period)

mortality$status <- NA
for (i in 1:nrow(mortality)) {
  if (mortality$mort_period[i] == "NA-NA") {
    mortality$status[i] <- "a"
  } else if (mortality$mort_period[i] == mortality$period[i]) {
    mortality$status[i] <- "d"
  } else if (mortality$status[i - 1] %in% c("d", "dp") 
             & mortality$tree.[i] == mortality$tree.[i - 1]) {
    mortality$status[i] <- "dp"
  } else {
    mortality$status[i] <- "a"
  }
}


# -- MORTALITY RATES: ALL SPECIES POOLED -----------------------------------

topocolor <- extra_plot_info %>%
  dplyr::rename(plot_group = plot) %>%
  dplyr::mutate(color = case_when(
    topo == "Xeric" ~ "brown",
    topo == "Xeric/Mesic" ~ "olivedrab4",
    topo == "Mesic" ~ "green",
    topo == "Hydric" ~ "blue"
  )) %>%
  dplyr::arrange(topo)

# Order moisture classes wet to dry, for consistent plot ordering below.
topocolor <- dplyr::bind_rows(
  topocolor %>% dplyr::filter(topo == "Hydric"),
  topocolor %>% dplyr::filter(topo == "Mesic"),
  topocolor %>% dplyr::filter(topo == "Xeric/Mesic"),
  topocolor %>% dplyr::filter(topo == "Xeric")
)

mortrates <- mortality %>%
  dplyr::group_by(plot_group, period, status) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  tidyr::spread(., status, n) %>%
  dplyr::filter(., !period %in% c("1982-1982", "1983-1983", "1986-1986", 
                                  "2016-2016", "NA-NA")) %>%
  dplyr::left_join(., PP_plot_data, by = c("plot_group" = "plot"))

# Install years for plot groups not resolved by the join above.
mortrates$install_year[mortrates$plot_group == "BL6_gap"] <- 1982
mortrates$install_year[mortrates$plot_group == "MRS4_gap"] <- 1983
mortrates$install_year[mortrates$plot_group == "MRS8-9-10"] <- 1986
mortrates$install_year[mortrates$plot_group == "MRS1_gap"] <- 1983

mortrates <- mortrates %>%
  dplyr::filter(!(install_year == 1986 & period == "1984-1986")) %>% # drop the census period before this plot group existed
  dplyr::left_join(., unique(PP_tree_data[, c("period", "duration")])) %>%
  dplyr::left_join(., extra_plot_info, by = c("plot_group" = "plot")) %>%
  dplyr::mutate(perc_mort = d / a) %>%
  dplyr::mutate(perc_mort_ann = perc_mort / duration)

mortrates$plot_group <- factor(mortrates$plot_group, 
                               levels = topocolor$plot_group)
mortrates$PercAnnMort <- mortrates$perc_mort_ann * 100

# Long (multi-decade) census intervals aren't comparable to the shorter
# ones, so exclude them here.
mortrates_sub <- subset(mortrates, duration < 4)

meanmort <- summarySE(as.data.frame(ungroup(mortrates_sub)), 
                      measurevar = "PercAnnMort", 
                      groupvars = c("period", "topo"))


# Quick look at mortality by topographic moisture class, all species
# pooled, before moving to the species-level breakdown used in the
# dashboard figures.
ggplot(meanmort, aes(x = period, y = PercAnnMort, colour = topo)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = PercAnnMort - se, ymax = PercAnnMort + se), 
                width = 0, alpha = 0.5) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "period", y = "% Annual Tree Mortality", colour = "Site Type") +
  theme(
    legend.position = c(0.14, 0.85),
    legend.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "line"),
    axis.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18)
  )


# -- MORTALITY RATES: BY SPECIES -------------------------------------------
# Same approach as above, but broken out by species and corrected to use
# dead / (dead + alive) as the mortality denominator, rather than
# dead / alive.

mortrates <- mortality %>%
  dplyr::group_by(plot_group, period, status, spec) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  tidyr::spread(., status, n) %>%
  dplyr::filter(., !period %in% c("1982-1982", "1983-1983", "1986-1986", 
                                  "2016-2016", "NA-NA")) %>%
  dplyr::left_join(., PP_plot_data, by = c("plot_group" = "plot"))

# If a plot/period/species has alive (a) trees but no dead (d) trees
# recorded, that means zero died -- fill in explicit 0s rather than NA.
mortrates <- mortrates %>%
  dplyr::mutate(d = ifelse(is.na(d) & !is.na(a), 0, d))

mortrates$install_year[mortrates$plot_group == "BL6_gap"] <- 1982
mortrates$install_year[mortrates$plot_group == "MRS4_gap"] <- 1983
mortrates$install_year[mortrates$plot_group == "MRS8-9-10"] <- 1986
mortrates$install_year[mortrates$plot_group == "MRS1_gap"] <- 1983

mortrates <- mortrates %>%
  dplyr::filter(!(install_year == 1986 & period == "1984-1986")) %>%
  dplyr::left_join(., unique(PP_tree_data[, c("period", "duration")])) %>%
  dplyr::left_join(., extra_plot_info, by = c("plot_group" = "plot")) %>%
  dplyr::mutate(perc_mort = d / (d + a)) %>%
  dplyr::mutate(perc_mort_ann = perc_mort / duration)

mortrates$plot_group <- factor(mortrates$plot_group, 
                               levels = topocolor$plot_group)
mortrates$PercAnnMort <- mortrates$perc_mort_ann * 100
mortrates_sub <- subset(mortrates, duration < 4)

# A handful of plots have PercAnnMort == NA because nothing was alive at
# the start of the period (only dead trees recorded) -- these are the
# "single last tree died" cases and are treated as 100% mortality below.
mortrates_all <- mortrates %>%
  mutate(PercAnnMort_adj = ifelse((is.na(a) & !is.na(d)), 1, PercAnnMort))

meanmortsp <- summarySE(as.data.frame(ungroup(mortrates_sub)), 
                        measurevar = "PercAnnMort", 
                        groupvars = c("period", "topo", "spec"))


meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
meanmortsp$spec <- as.factor(meanmortsp$spec)
meanmortsp <- subset(meanmortsp, spec != "PIFL")  # too few Limber Pine records to summarize reliably

facet_names <- c(
  `ABLA` = "Subalpine Fir", `PICO` = "Lodgepole Pine",
  `PIEN` = "Engelmann Spruce", `POTR` = "Quaking Aspen"
)

# Quick look at mortality by species, all plots, before restricting to
# the species/plot subset used in the dashboard figures.
ggplot(meanmortsp, aes(x = period, y = PercAnnMort, colour = topo)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = PercAnnMort - se, ymax = PercAnnMort + se), 
                width = 0, alpha = 0.5) +
  labs(x = "period", y = "% Annual Tree Mortality", colour = "Site Type") +
  facet_wrap(~spec, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme(
    legend.position = "right",
    legend.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "line"),
    axis.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18)
  )


# -- RESTRICT TO RECOMMENDED SPECIES/PLOTS -----------------------------------
# Per R. Andrus: too few Quaking Aspen (POTR) trees for confidence, and
# Lodgepole Pine (PICO) mortality is largely density-dependent rather than
# climate-driven, so both are dropped from the dashboard figures. Also
# restricted to the plots where climate (rather than density dependence)
# is thought most likely to drive mortality patterns.

mortrates_sub <- subset(mortrates_sub, plot_group %in% 
                          c("BW3", "MRS4", "MRS5", "BL6", "MRS7"))

mortrates_sub_all <- mortrates_sub %>%
  mutate(PercAnnMort_adj = ifelse((is.na(a) & !is.na(d)), 1, PercAnnMort))

meanmortsp <- summarySE(
  as.data.frame(ungroup(mortrates_sub %>% 
                          filter(plot_group %in% c("BW3", "MRS4", "MRS5", "BL6",
                                                   "MRS7")))),
  measurevar = "PercAnnMort", groupvars = c("period", "topo", "spec")
)

meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
meanmortsp$spec <- as.factor(meanmortsp$spec)
meanmortsp <- subset(meanmortsp, spec == "ABLA" | spec == "PIEN")


# ============================================================================
# FIGURE 1: % annual mortality by census period and topographic class
# ============================================================================
# Subalpine Fir (ABLA) and Engelmann Spruce (PIEN) only, restricted to the
# plots recommended by R. Andrus. Points show the mean ± SE across plots
# within a topographic moisture class; the linear trend line is fit to
# Subalpine Fir only.

g1 <- ggplot(meanmortsp, aes(x = period, y = PercAnnMort, colour = spec)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = PercAnnMort - se, ymax = PercAnnMort + se), 
                width = 0, alpha = 0.5) +
  labs(x = "Census period", y = "Annual Tree Mortality (%)", colour = NULL) +
  scale_colour_discrete(labels = c("Subalpine Fir", "Engelmann Spruce")) +
  geom_smooth(method = lm, se = FALSE, 
              data = subset(meanmortsp, spec == "ABLA")) +
  theme(
    legend.position = c(0.20, 0.75),
    legend.background = element_rect(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(2, "line"),
    axis.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = -90, hjust = 0)
  ) +
  facet_wrap(~topo)

ggsave(g1, file = file.path(plot_dir, "subalpine_tree_mortality.jpg"), 
       height = 5, width = 12)


# ============================================================================
# FIGURE 2: % annual mortality anomaly by species and topographic moisture class
# ============================================================================
# Companion figure to Figure 1: rather than points with error bars, this
# shows each period's mortality rate as a bar, colored by whether it's
# above or below that species/moisture class long-term mean --
# Mortality is weighted by the number of trees alive+dead in that plot/period, 
# so plots with more trees contribute proportionally more to the average.

meanmortsp_adj <- mortrates_sub_all %>%
  filter(plot_group %in% c("BW3", "MRS4", "MRS5", "BL6", "MRS7")) %>%
  mutate(
    ct_a = ifelse(!is.na(a), a, 0),
    ct_d = ifelse(!is.na(d), d, 0),
    sum_tree = ct_a + ct_d
  ) %>%
  ungroup() %>%
  filter(sum_tree != 0) %>%
  dplyr::group_by(topo, spec, period) %>%
  dplyr::summarize(PercAnnMort = weighted.mean(PercAnnMort_adj, sum_tree),
                   .groups = "drop") %>%
  tidyr::separate(., col = "period", into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start), end = as.numeric(end), 
         year = (start + end) / 2) %>%
  filter(spec %in% c("ABLA", "PIEN")) %>%
  dplyr::group_by(spec, topo) %>%
  dplyr::mutate(mean_mort = mean(PercAnnMort, na.rm = TRUE), 
                .groups = "drop") %>%
  mutate(pos_neg = factor(ifelse(PercAnnMort > mean_mort, "pos", "neg"), 
                          levels = c("neg", "pos")))


facet_names <- c(
  `ABLA` = "Subalpine Fir", `PICO` = "Lodgepole Pine",
  `PIEN` = "Engelmann Spruce", `POTR` = "Quaking Aspen",
  `Hydric` = "Hydric", `Mesic` = "Mesic", `Xeric` = "Xeric"
)

g2 <- ggplot(
  meanmortsp_adj %>% filter(spec %in% c("ABLA", "PIEN")),
  aes(x = year, y = PercAnnMort, fill = pos_neg)
) +
  geom_col() +
  labs(x = "Census period", y = "Annual Tree Mortality (%)", colour = NULL) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  facet_grid(spec ~ topo, labeller = as_labeller(facet_names)) +
  theme(legend.position = "none")

ggsave(g2,
       file = file.path(plot_dir, "subalpine_tree_mortality_long-term_mean.jpg"),
       height = 8, width = 10, scale = 0.5
)

