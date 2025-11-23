# Script to download and plot c1/d1 air temp and ppt trends by quarter
# and/or by winter/summer
# and/or annual
# Updated 21 Nov 2025

# -- SETUP ------
# Clean environment
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

options(stringsAsFactors = F)

na_vals <- c("NP", "NA", NA, "NaN", NaN, ".")

library(tidyverse)
library(EDIutils) # Handy tools for interacting with EDI's API
library(trend) # for Mann-Kendall and Theil-Sen tests
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(lemon)
library(ggthemes)
theme_set(theme_bw())
library(viridisLite) #for color-blind accessibility

# only need to download once
download_data <- TRUE

# download data -----------------------------------------------------------
# Note: if you have already downloaded SOME data, the read_data_package_archive
# function will not work as it does not want to overwrite. Clear your /data
# directories and then return to run the following code.
if (download_data) {
  # download the data from EDI
  # 184 is c1 ppt
  # 186 is d1 ppt
  # 315 is the sdl ppt
  # 314 is sdl temp
  # 185 is c1 temp
  # 187 is d1 temp
  
  data_dir <- file.path("c1_d1_sdl_temp_ppt", "data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  scope <- "knb-lter-nwt" # Niwot scope

  # Note: the overwrite argument does not work, so clear out any existing
  # copies before running this
  for (id in c(
    "184", "186", "315",
    "314", "185", "187"
  )) {
    # Ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")

    # Display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")

    # Download the data
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }

  # Update the below so you remember to cite the data correctly
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for C1 chart recorder, 1952 - ongoing. ver 8. Environmental Data Initiative. https://doi.org/10.6073/pasta/5f896d7d1eb19649ee2ff9fca7166160. Accessed 2025-07-09."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for D1 chart recorder, 1952 - ongoing. ver 6. Environmental Data Initiative. https://doi.org/10.6073/pasta/02b125f7d0dc87a6d40e6ca7001471ce. Accessed 2025-07-09."
  # [1] "White, C., J. Morse, K. Chowanski, T. Kittel, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for Saddle, 1981 - ongoing. ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/8caa444276521b797e4a16485dc0e137. Accessed 2025-07-11."
  # [1] "White, C., J. Morse, H. Brandes, K. Chowanski, T. Kittel, M. Losleben, and M. Moore. 2025. Homogenized, gap-filled, daily air temperature data for Saddle, 1986 - ongoing. ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/49a2171ec028c4ddef298752bc9ebe8d. Accessed 2025-07-11."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for C1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/0b4f8747ea72258e46681bb511c262f3. Accessed 2025-07-11."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for D1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/b05689181c21ac40d35b6c1c01e2f8a5. Accessed 2025-07-11."
  # Unzip all zip files in the directory
  # This overwrites the manifests, but you don't really need them.
  for (fname in list.files(data_dir,
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}

# functions and plot formats-----------------------------------------------------------

# define function to perform
perform_trend_analysis <- function(data, value_col = "tot_ppt", year_col = "year") {
  # Remove any missing values
  clean_data <- data[!is.na(data[[value_col]]), ]
  
  if (nrow(clean_data) < 3) {
    return(data.frame(
      mk_tau = NA,
      mk_pvalue = NA,
      mk_trend = "insufficient data",
      ts_slope = NA,
      ts_intercept = NA
    ))
  }
  # Ensure data is sorted by year before trend analysis
  data <- data %>% arrange(!!sym(year_col))
  # Mann-Kendall test
  mk_test <- mk.test(data[[value_col]])
  
  # Theil-Sen slope
  ts_test <- as.numeric(sens.slope(data[[value_col]])$estimates)
  
  # Calculate intercept
  intercepts <- data[[value_col]] - ts_test * data[[year_col]]
  sens_intercept <- median(intercepts)
  
  # Determine trend direction
  trend_direction <- case_when(
    mk_test$p.value > 0.05 ~ "no trend",
    mk_test$statistic > 0 ~ "increasing",
    mk_test$statistic < 0 ~ "decreasing",
    TRUE ~ "no trend"
  )
  
  return(data.frame(
    mk_tau = as.numeric(mk_test$statistic),
    mk_pvalue = as.numeric(mk_test$p.value),
    mk_trend = trend_direction,
    ts_slope = ts_test,
    ts_intercept = sens_intercept
  ))
}

# Create plot for a specific season with all sites
create_seasonal_plot <- function(seasonal_data, trend_results, target_season,
                                 y_col = "tot_ppt", y_label = "Precipitation (mm)",
                                 plot_title = NULL, slope_units = "mm/yr") {
  # Filter data for the target season
  plot_data <- seasonal_data %>%
    filter(season == target_season) %>%
    left_join(trend_results %>% filter(season == target_season),
              by = c("site", "season")
    ) %>%
    mutate(
      # Create significance flag for line types
      significant = ifelse(mk_pvalue < 0.05, "Significant", "Non-significant"),
      # Handle cases where significance is NA
      significant = ifelse(is.na(significant), "Insufficient data", significant)
    )
  
  # Define colors for sites (you may need to adjust these based on your sites)
  site_colors <- c(
    "C1" = "#D55E00",
    "D1" = "#56B4E9",
    "SDL" = "#009E73"
  )
  site_linetypes <- c(
    "C1" = "solid",
    "D1" = "dashed",
    "SDL" = "dotted"
  )
  
  # Calculate trend lines for each site's data range
  trend_lines <- plot_data %>%
    group_by(site) %>%
    filter(!is.na(ts_slope) & !is.na(ts_intercept)) %>%
    summarize(
      min_year = min(year, na.rm = TRUE),
      max_year = max(year, na.rm = TRUE),
      ts_slope = first(ts_slope),
      ts_intercept = first(ts_intercept),
      significant = first(significant),
      .groups = "drop"
    ) %>%
    mutate(
      y_start = ts_intercept + ts_slope * min_year,
      y_end = ts_intercept + ts_slope * max_year
    )%>%
    filter(significant=="Significant") #this will keep only significant trends; a line will be drawn for only significant trends
  
  # Create the plot title
  if (is.null(plot_title)) {
    plot_title <- paste(stringr::str_to_title(target_season), "Trends by Site")
  }
  
  # Create main plot using tidy evaluation for dynamic y column
  p <- ggplot(plot_data, aes(x = year, y = !!sym(y_col), color = site)) +
    # geom_point(alpha = 0.6, size = 1.5) +
    geom_line(alpha = 0.7, size = 1) +
    geom_segment(
      data = trend_lines,
      aes(
        x = min_year, y = y_start,
        xend = max_year, yend = y_end,
        color = site, linetype = site
      ),
      size = 1.2
    ) +
    scale_color_manual(values = site_colors, name = "Site") +
    scale_linetype_manual(values = site_linetypes)+
    #scale_linetype_manual(
    # values = c(
    #  "Significant" = "solid",
    #  "Non-significant" = "dashed",
    # "Insufficient data" = "dotted"
    # ),
    # name = "Mann-Kendall\nSignificance"
    # ) +
    ggtitle(plot_title) +
    ylab(y_label) +
    xlab("Year") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(linetype = "solid")))
  
  # Create table with slope information
  slope_data <- plot_data %>%
    select(site, ts_slope, mk_pvalue, mk_trend) %>%
    distinct() %>%
    mutate(
      slope_round = round(ts_slope, 3),
      p_round = round(mk_pvalue, 3),
      significance = ifelse(mk_pvalue < 0.05, "*", ""),
      display_text = paste0(slope_round, significance)
    )
  
  # Create Grid Graphical Object (grob) with dynamic slope units
  inset_data <- data.frame(
    Site = slope_data$site,
    stringsAsFactors = FALSE
  )
  
  #Rename sites for display
  inset_data <- inset_data %>%
    mutate(Site = case_when(
      Site == "C1" ~ "Subalpine (C1)",
      Site == "D1" ~ "Alpine (D1)",
      Site == "SDL" ~ "Saddle (SDL)",  # Keep SDL as is, or change to something else
      TRUE ~ Site  # Keep any other sites as is
    ))
  
  inset_data[[paste0("Slope (", slope_units, ")")]] <- slope_data$display_text
  inset_data[["p-value"]] <- slope_data$p_round
  
  # Create table grob
  table_grob <- tableGrob(inset_data,
                          rows = NULL,
                          theme = ttheme_minimal(
                            core = list(
                              fg_params = list(cex = 0.8),
                              bg_params = list(fill = "white", alpha = 0.8)
                            ),
                            colhead = list(
                              fg_params = list(cex = 0.9, fontface = "bold"),
                              bg_params = list(fill = "lightgray", alpha = 0.8)
                            )
                          )
  )
  
  # Add colored rectangles for sites in table
  for (i in 1:nrow(inset_data)) {
    site_name <- inset_data$Site[i]
    #Map display names back to original site codes for color lookup
    original_site <- case_when(
      site_name == "Subalpine (C1)" ~ "C1",
      site_name == "Alpine (D1)" ~ "D1",
      site_name =="Saddle (SDL)" ~"SDL",
      TRUE ~ site_name
    )
    color <- site_colors[original_site]
    table_grob <- gtable_add_grob(table_grob,
                                  rectGrob(gp = gpar(fill = color, alpha = 0.3, col = NA)),
                                  t = i + 1, l = 1, r = 1
    )
  }
  
  #Add the table with control over spacing
  p_with_table <- grid.arrange(
    p,
    table_grob,
    nrow = 2,
    heights = c(5, 1)  # Ratio of 3:1 (plot taller than table)
  )
  
  return(p_with_table)
  
  # Add the table as annotation - use dynamic y column for positioning
  #y_values <- plot_data[[y_col]]
  #p_with_inset <- p +
  # annotation_custom(table_grob,
  #  xmin = max(plot_data$year) - (max(plot_data$year) - min(plot_data$year)) * 0.35,
  # xmax = max(plot_data$year) - (max(plot_data$year) - min(plot_data$year)) * 0.05,
  #ymin = max(y_values, na.rm = TRUE) - (max(y_values, na.rm = TRUE) - min(y_values, na.rm = TRUE)) * 0.4,
  #ymax = max(y_values, na.rm = TRUE) - (max(y_values, na.rm = TRUE) - min(y_values, na.rm = TRUE)) * 0.05
  #)
  
  #return(p_with_inset)
}


# Read and merge ppt and temp datasets --------------------------------------------------

ppt_seasonal <- read_csv(file.path( data_dir, "sdl_daily_precip_gapfilled.cw.data.csv")) %>%
  bind_rows(., read_csv(file.path( data_dir, "d1_infilled_precip_daily.tk.data.csv"))) %>%
  bind_rows(., read_csv(file.path( data_dir, "c1_infilled_precip_daily.tk.data.csv"))) %>%
  mutate(
    month = lubridate::month(date),
    year = lubridate::year(date)
  ) %>%
  mutate(
    season =
      case_when(
        month %in% c(6, 7, 8) ~ "summer",
        TRUE ~ "winter"
      )
  ) %>%
  rename(site = local_site) %>%
  group_by(site, year, season) %>%
  summarize(tot_ppt = sum(precip), .groups = "drop") %>%
  filter(!(site == "SDL")) # %>%
# We have dropped sdl entirely from ppt because it's not reliable enough even
# after this:
# filter(!(site == "SDL" & year < 1990))

temp_seasonal <- read_csv(file.path( data_dir, "d1_infilled_temp_daily.tk.data.csv"),
  guess_max = 100000
) %>%
  bind_rows(., read_csv(file.path( data_dir, "c1_infilled_temp_daily.tk.data.csv"),
    guess_max = 100000
  )) %>%
  bind_rows(., read_csv(file.path( data_dir, "sdl_daily_airtemp_gapfilled.cw.data.csv"),
    guess_max = 100000
  ) %>%
    rename(
      min_temp = airtemp_min_homogenized,
      max_temp = airtemp_max_homogenized
    )) %>%
  mutate(
    month = lubridate::month(date),
    year = lubridate::year(date)
  ) %>%
  mutate(
    season =
      case_when(
        month %in% c(6, 7, 8) ~ "summer",
        month %in% c(9, 10, 11) ~ "fall",
        month %in% c(12, 1, 2) ~ "winter",
        month %in% c(3, 4, 5) ~ "spring",
        TRUE ~ NA
      ),
    mean_temp = (min_temp + max_temp) / 2
  ) %>%
  rename(site = local_site) %>%
  group_by(site, year, season) %>%
  summarize(mean_temp = mean(mean_temp, na.rm = FALSE), .groups = "drop")

# Create annual aggregations for D1 site ---------------------------------------
ppt_annual_d1 <- ppt_seasonal %>%
  filter(site == "D1") %>%
  group_by(site, year) %>%
  summarize(tot_ppt = sum(tot_ppt, na.rm = TRUE), .groups = "drop")

# Recalculate annual temperature from daily data to account for different season lengths
temp_annual_d1 <- read_csv(file.path( file.path( data_dir, "d1_infilled_temp_daily.tk.data.csv")),
  guess_max = 100000
) %>%
  mutate(
    month = lubridate::month(date),
    year = lubridate::year(date),
    mean_temp = (min_temp + max_temp) / 2
  ) %>%
  rename(site = local_site) %>%
  filter(site == "D1") %>%
  group_by(site, year) %>%
  summarize(mean_temp = mean(mean_temp, na.rm = TRUE), .groups = "drop")

# Calculate trends by site and season-------------------------------------------
trend_results_ppt <- ppt_seasonal %>%
  # drop sdl not reliable
  # filter(site!= 'SDL') %>%
  group_by(site, season) %>%
  do(perform_trend_analysis(., value_col = "tot_ppt", year_col = "year")) %>%
  ungroup()

trend_results_temp <- temp_seasonal %>%
  group_by(site, season) %>%
  do(perform_trend_analysis(., value_col = "mean_temp", year_col = "year")) %>%
  ungroup()

# Calculate annual trends for D1 site ------------------------------------------
trend_results_ppt_annual_d1 <- ppt_annual_d1 %>%
  do(perform_trend_analysis(., value_col = "tot_ppt", year_col = "year")) %>%
  mutate(site = "D1", season = "annual")

trend_results_temp_annual_d1 <- temp_annual_d1 %>%
  do(perform_trend_analysis(., value_col = "mean_temp", year_col = "year")) %>%
  mutate(site = "D1", season = "annual")

# Calculate trends by site and season-------------------------------------------
# Generate individual plots
summer_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "summer",
  plot_title = "Jun-Aug"
)
winter_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "winter",
  plot_title = "Sep-May"
)
summer_plot_ppt

summer_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "summer",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = "Summer\n(Jun-Aug)", slope_units = "degC/yr"
)
summer_plot_temp

spring_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "spring",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = "Spring\n(Mar-May)", slope_units = "degC/yr"
)
spring_plot_temp

fall_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "fall",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = "Fall\n(Sep-Nov)", slope_units = "degC/yr"
)
fall_plot_temp

winter_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "winter",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = "Winter\n(Dec-Feb)", slope_units = "degC/yr"
)
winter_plot_temp

summer_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "summer",
  plot_title = "Jun-Aug"
)
summer_plot_ppt

winter_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "winter",
  plot_title = "Sep-May"
)
winter_plot_ppt

# Combined plots ---------------------------------------------------------------

# Extract legend from one plot
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

shared_legend_temp <- get_legend(spring_plot_temp)
shared_legend_ppt <- get_legend(summer_plot_ppt)


jpeg("c1_d1_sdl_temp_ppt/figures/combined_plot_temp.jpg", width = 12, height = 10, units = "in", res = 300)
grid.arrange(winter_plot_temp + theme(legend.position = "none"),
  spring_plot_temp + theme(legend.position = "none"),
  summer_plot_temp + theme(legend.position = "none"),
  fall_plot_temp + theme(legend.position = "none"),
  shared_legend_temp,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
  heights = c(1, 1, 0.1)
)
dev.off()

jpeg("c1_d1_sdl_temp_ppt/figures/combined_plot_ppt.jpg", width = 12, height = 6, units = "in", res = 300)
combined_plot_ppt <- grid.arrange(winter_plot_ppt + theme(legend.position = "none"),
  summer_plot_ppt + theme(legend.position = "none"),
  shared_legend_ppt,
  ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1, 2), c(5, 5)),
  heights = c(1, 0.1)
)
dev.off()

# Summer temp simple plot ------------------------------------------------------
# create the thiel sen slope segments
trend_segments <- trend_results_temp %>%
  mutate(
    # Create significance flag for line types
    significant = ifelse(mk_pvalue < 0.05, "Significant", "Non-significant"),
    # Handle cases where significance is NA
    significant = ifelse(is.na(significant), "Insufficient data", significant)
  ) %>%
  group_by(site, season) %>%
  filter(!is.na(ts_slope) & !is.na(ts_intercept)) %>%
  inner_join(temp_seasonal %>% filter(season == "summer"),
    by = c("site", "season")
  ) %>%
  summarize(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    ts_slope = first(ts_slope),
    ts_intercept = first(ts_intercept),
    significant = first(significant),
    .groups = "drop"
  ) %>%
  mutate(
    y_start = ts_intercept + ts_slope * min_year,
    y_end = ts_intercept + ts_slope * max_year
  )

# simpler plots


###
slope_data <- trend_results_temp %>%
  filter(site %in% c("C1", "D1") & season %in% c("summer"))

# Create labels with slope values
c1_slope <- slope_data$ts_slope[slope_data$site == "C1"]
d1_slope <- slope_data$ts_slope[slope_data$site == "D1"]

legend_labels <- c(
  "C1" = paste0("Subalpine +", round(c1_slope, 3), " °C per year"),
  "D1" = paste0("Tundra + ", round(d1_slope, 3), " °C per year")
)


c1_d1_temps <- ggplot(
  temp_seasonal %>% filter(site %in% c("C1", "D1") &
    season %in% c("summer")),
  aes(x = year, y = mean_temp, color = site)
) +
  geom_point() +
  geom_line() +
  geom_segment(
    data = trend_segments %>%
      filter(site %in% c("C1", "D1") &
        season %in% c("summer") &
        significant == "Significant"),
    aes(
      x = min_year, y = y_start,
      xend = max_year, yend = y_end,
      color = site
    ),
    size = 1.2
  ) +
  scale_color_manual(
    values = c("C1" = "red", "D1" = "blue"),
    labels = legend_labels,
    name = NULL
  ) +
  labs(
    title = "Summer mean temperature on Niwot Ridge by biome",
    x = "Year",
    y = "Temperature (°C)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.margin = margin(6, 6, 6, 6)
  )

ggsave(c1_d1_temps,
  filename = "c1_d1_sdl_temp_ppt/figures/c1_d1_summer_mean_temp.jpg",
  width = 12, height = 8, units = "in", dpi = 300,
  scale = 0.5
)

# Anomaly plots ----------------------------------------------------------------

# Seasonal Anomalies for D1 Precipitation and Temperature------------
# Create seasonal anomaly data for temperature
anom_temp_season <- temp_seasonal %>%
  group_by(site, season) %>%
  mutate(
    mean_temp_season = mean(mean_temp, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(
    anom_temp = mean_temp - mean_temp_season,
    posneg = ifelse(anom_temp > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

# Create seasonal anomaly data for precipitation
anom_ppt_season <- ppt_seasonal %>%
  group_by(site, season) %>%
  mutate(
    mean_ppt_season = mean(tot_ppt, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(
    # anom_ppt = tot_ppt/mean_ppt_season,
    anom_ppt = (tot_ppt - mean_ppt_season) / mean_ppt_season * 100,
    posneg = ifelse(anom_ppt > 1, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

#Create seasonal anomaly plot for temperature
anom_temp_D1 <- ggplot(anom_temp_season %>%
                         filter(site == "D1") %>%
                         ungroup() %>%
                         mutate(
                           season =
                             fct_relevel(
                               season, "spring", "summer",
                               "fall", "winter"
                             )
                         ), aes(x = year, y = anom_temp)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("pos" = "#99000d", "neg" = "#034e7b")) +
  labs(y = "Seasonal temperature anomaly (°C)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = identity, breaks = NULL, name = expression(hotter %<->% colder))) +
  theme_hc() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  facet_wrap(~season, labeller = labeller(season = c("spring" = "Spring (Mar-May)", "summer" = "Summer (Jun-Aug)", "fall" = "Fall (Sep-Nov)", "winter" = "Winter (Dec-Feb)")))

ggsave(anom_temp_D1,
       filename = "c1_d1_sdl_temp_ppt/figures/d1_temp_anom_by_season.png",
       scale = 0.5, width = 8, height = 8
)

#Create seasonal anomaly plot for precipitation
anom_ppt_D1 <- ggplot(anom_ppt_season %>%
  ungroup() %>%
  filter(site == "D1") %>%
  mutate(
    season =
      fct_relevel(
        season, "summer",
        "winter"
      )
  ), aes(x = year, y = anom_ppt)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("pos" = "#034e7b", "neg" = "#99000d")) +
  labs(y = "Seasonal precipitation anomaly (%)", x = "") +
  scale_y_symmetric(sec.axis = sec_axis(trans = identity, breaks = NULL, name = expression(wetter %<->% drier))) +
  theme_hc() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  ) +
  facet_wrap(~season, labeller = labeller(season = c("summer" = "Jun-Aug", "winter" = "Sep-May")))

ggsave(anom_ppt_D1,
  file = "c1_d1_sdl_temp_ppt/figures/d1_ppt_anom_by_season.png",
  scale = 0.5, width = 8, height = 6
)


# Annual anomaly plots for D1 precipitation and temperature ---------------------------------------------

# Create annual anomaly data for temperature
anom_temp_annual_d1 <- temp_annual_d1 %>%
  mutate(
    mean_temp_overall = mean(mean_temp, na.rm = TRUE),
    anom_temp = mean_temp - mean_temp_overall,
    posneg = ifelse(anom_temp > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

# Create annual anomaly data for precipitation
anom_ppt_annual_d1 <- ppt_annual_d1 %>%
  mutate(
    mean_ppt_overall = mean(tot_ppt, na.rm = TRUE),
    anom_ppt = (tot_ppt - mean_ppt_overall) / mean_ppt_overall * 100,
    posneg = ifelse(anom_ppt > 0, "pos", "neg") %>%
      factor(c("pos", "neg"))
  )

# Create annual temperature anomaly plot for D1
anom_temp_annual_d1_plot <- ggplot(anom_temp_annual_d1, aes(x = year, y = anom_temp)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("pos" = "#99000d", "neg" = "#034e7b")) +
  labs(
    # title = "Annual Temperature Anomalies - D1 Site",
    y = "Annual temperature anomaly (°C)",
    x = "Year"
  ) +
  scale_y_symmetric(sec.axis = sec_axis(trans = identity, breaks = NULL, name = expression(hotter %<->% colder))) +
  theme_hc() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Create annual precipitation anomaly plot for D1
anom_ppt_annual_d1_plot <- ggplot(anom_ppt_annual_d1, aes(x = year, y = anom_ppt)) +
  geom_col(aes(fill = posneg)) +
  scale_fill_manual(values = c("pos" = "#034e7b", "neg" = "#99000d")) +
  labs(
    y = "Annual precipitation anomaly (%)",
    x = "Year"
  ) +
  scale_y_symmetric(sec.axis = sec_axis(
    trans = "identity", breaks = NULL,
    name = expression(wetter %<->% drier)
  )) +
  theme_hc() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Save annual anomaly plots
ggsave(anom_temp_annual_d1_plot,
  filename = "c1_d1_sdl_temp_ppt/figures/d1_temp_anom_annual.png",
  scale = 0.5, width = 8, height = 6
)

ggsave(anom_ppt_annual_d1_plot,
  filename = "c1_d1_sdl_temp_ppt/figures/d1_ppt_anom_annual.png",
  scale = 0.5, width = 8, height = 6
)
