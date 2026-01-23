# Analyzes and visualizes subalpine data from the University Camp snotel site @ 10,330ft. Figures visualize SWE (max day and mm),  day of last snowmelt, and fraction of melt pre max SWE.
# Last updated 1/21/2026

# setup
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(snotelr)
library(EDIutils) # Handy tools for interacting with EDI's API
library(trend) # for Mann-Kendall and Theil-Sen tests
library(gridExtra)
library(gtable)
library(grid)


# only need to download once
download_data_snotel <- TRUE
download_data_c1 <- FALSE

data_dir <- file.path("c1_snow", "data")
figures_dir <- file.path("c1_snow", "figures")

# Create directories if they don't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# functions -----------------------------------------------------------
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
  data <- data %>% arrange(year_col)
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
                                 plot_title = NULL, slope_units = "mm/yr",
                                 include_table = TRUE) {
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
      significant = ifelse(is.na(significant), "Insufficient data", significant),
      # Make it a factor with all levels defined
      significant = factor(significant, levels = c("Significant", "Non-significant", "Insufficient data"))
    )

  # Define colors for sites - using consistent color scheme
  site_colors <- c(
    "c1" = "#D55E00"
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
      y_end = ts_intercept + ts_slope * max_year,
      # Make it a factor with all levels defined
      significant = factor(significant, levels = c("Significant", "Non-significant", "Insufficient data"))
    )

  # Create dummy data to ensure both significance levels appear in legend
  dummy_year <- min(plot_data$year, na.rm = TRUE)
  dummy_y <- min(plot_data[[y_col]], na.rm = TRUE)

  dummy_trend <- data.frame(
    min_year = rep(dummy_year, 3),
    max_year = rep(dummy_year, 3), # Same as min_year so length = 0
    y_start = rep(dummy_y, 3),
    y_end = rep(dummy_y, 3), # Same as y_start so length = 0
    ts_slope = rep(0, 3),
    ts_intercept = rep(0, 3),
    site = rep(plot_data$site[1], 3),
    significant = factor(c("Significant", "Non-significant", "Insufficient data"),
      levels = c("Significant", "Non-significant", "Insufficient data")
    )
  )

  # Add dummy data to trend_lines
  trend_lines <- bind_rows(trend_lines, dummy_trend)

  # Create the plot title
  if (is.null(plot_title)) {
    plot_title <- paste(stringr::str_to_title(target_season), "Trends by Site")
  }

  # Create main plot using tidy evaluation for dynamic y column
  p <- ggplot(plot_data, aes(x = year, y = !!sym(y_col), color = site)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_line(linewidth = 1) +
    geom_segment(
      data = trend_lines,
      aes(
        x = min_year, y = y_start,
        xend = max_year, yend = y_end,
        color = site, linetype = significant
      ),
      linewidth = 1.2
    ) +
    scale_color_manual(
      values = site_colors, name = "Site",
      labels = c(
        "c1" = "Subalpine"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Significant" = "solid",
        "Non-significant" = "dashed",
        "Insufficient data" = "dotted"
      ),
      name = "Mann-Kendall\nSignificance",
      limits = c("Significant", "Non-significant", "Insufficient data"),
      drop = FALSE,
      breaks = c("Significant", "Non-significant", "Insufficient data")
    ) +
    ggtitle(plot_title) +
    ylab(y_label) +
    xlab("Year") +
    theme_classic() +
    theme(
      legend.text = element_text(size = 8),
      legend.position = if (!include_table) "none" else "bottom",
      legend.box = "horizontal",
      legend.box.background = element_rect(color = "black", linewidth = 0.5),
      legend.box.margin = margin(6, 6, 6, 6),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(
        order = 1, nrow = 1, override.aes = list(size = 1.2), keywidth = 1,
        title.theme = element_text(size = 7),
        label.theme = element_text(size = 7)
      ),
      linetype = guide_legend(
        order = 2, nrow = 3, override.aes = list(size = 1.2, linetype = c("solid", "dashed", "dotted")), keywidth = 2.5,
        title.theme = element_text(size = 7),
        label.theme = element_text(size = 7)
      )
    )

  # Return early if table is not needed
  if (!include_table) {
    return(p)
  }

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

  # Rename sites for display
  inset_data <- inset_data %>%
    mutate(Site = case_when(
      Site == "c1" ~ "Subalpine",
      TRUE ~ Site # Keep any other sites as is
    ))

  inset_data[[paste0("Slope (", slope_units, ")")]] <- slope_data$display_text
  inset_data[["p-value"]] <- slope_data$p_round

  # Create table Grob
  table_grob <- tableGrob(inset_data,
    rows = NULL,
    theme = ttheme_default(
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
    # Map display names back to original site codes for color lookup
    original_site <- case_when(
      site_name == "Subalpine" ~ "c1",
      TRUE ~ site_name
    )
    color <- site_colors[original_site]
    table_grob <- gtable_add_grob(table_grob,
      list(rectGrob(gp = gpar(fill = NA, col = color, lwd = 3))),
      t = i + 1, l = 1, r = 1
    )
  }

  # Extract the legend from the plot
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend_grob <- tmp$grobs[[leg]]

  # Remove legend from the main plot
  p_no_legend <- p + theme(legend.position = "none")

  # Create spacers for centering
  spacer_left <- rectGrob(gp = gpar(col = NA, fill = NA))
  spacer_right <- rectGrob(gp = gpar(col = NA, fill = NA))

  # Create a combined grob with legend and table side by side, centered
  bottom_row <- arrangeGrob(spacer_left, legend_grob, table_grob, spacer_right,
    ncol = 4, widths = c(0.5, 3, 2, 0.5)
  )

  # Combine plot and bottom row
  p_with_table <- grid.arrange(
    p_no_legend,
    bottom_row,
    nrow = 2,
    heights = c(5, 1) # Ratio of 5:1 (plot taller than bottom elements)
  )
  return(p_with_table)
}

# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data_c1) {
  # download the data from EDI
  # 185 is c1 temp

  scope <- "knb-lter-nwt" # Niwot scope

  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c("185")) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")

    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")

    # download the data
    read_data_package_archive(packageID, path = data_dir)
    print(read_data_package_citation(packageID))
  }

  # update the below so you remember to cite it correctly
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for C1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/0b4f8747ea72258e46681bb511c262f3. Accessed 2025-07-11."

  # unzip all zip files in the directory
  # overwrites the manifests but don't really need them.
  for (fname in list.files(data_dir,
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = data_dir)
  }
}

if (download_data_snotel) {
  # site_meta_data <- snotel_info()
  # head(site_meta_data)

  # but we know it's 838, so feed in manually
  # snotel data?

  # meta_data <- snotel_info() #if you need to sort out
  # again other stations
  # data starts Oct 1 1978
  # but note there's no swe data til oct 1 1979

  # often you have to run this twice I am not sure why
  snotel_download(
    site_id = 838, path = data_dir,
    internal = F
  )
}

# read and munge snow  --------------------------------------------------
# Check if SNOTEL data file exists
snotel_file <- file.path(data_dir, "snotel_838.csv")
if (!file.exists(snotel_file)) {
  stop(paste0(
    "SNOTEL data file not found at: ", snotel_file, "\n",
    "Please either:\n",
    "  1. Set download_data_snotel <- TRUE at the top of the script to download the data, or\n",
    "  2. Make sure the 'snotel_838.csv' file is in the '", data_dir, "' directory"
  ))
}

snow_data_snotel <- read_csv(snotel_file)

# one missing value but it's july and 0 long before and after
snow_data_snotel <- snow_data_snotel %>%
  arrange(date) %>%
  # one missing value but it's july and 0 long before and after
  mutate(
    snow_water_equivalent = ifelse(date == "2022-07-16", 0, snow_water_equivalent)
  )


# check what dates avail
range(snow_data_snotel$date[!is.na(snow_data_snotel$snow_water_equivalent)]) # so only water year 1980 - 2025ish good for snow at this point


# calculate snow phenology,
# note this only runs if 365d in the year
# to the last/partial year will be cut off
snow_phenology_snotel <- snotel_phenology(snow_data_snotel)


# add the current year we are in if you wish
# max_2025<-max(snow_data_snotel$snow_water_equivalent[snow_data_snotel$date >= '2025-01-01'])
# max_2025_doy <- snow_data_snotel$date[snow_data_snotel$snow_water_equivalent == max_2025 &
#                         snow_data_snotel$date >= '2025-01-01'] %>%
#   lubridate::yday(.) %>% max(.)
#
# snow_phenology_snotel <- snow_phenology_snotel %>%
#   bind_rows(., data.frame (year = 2025, max_swe = max_2025, max_swe_doy = max_2025_doy))




# Calculate frac melt pre max SWE -------------------------------------------
# Rationale from NCC paper
# We introduce a metric derived from daily SWE observations that complements
# the date of maximum SWE and provides additional hydrologic information.
# The fraction of cumulative annual snow water resources that has melted
# before a given date i, FMi, was computed for each station-year. This daily
# metric was computed in three steps for each of two dates: 1) the date of
# maximum SWE and 2) April 1st. First, daily melt (Figure S1; blue bars)
# was computed as the daily decrease in SWE, presented as positive values.
# Second, cumulative daily melt (Figure S1; red hashed line) was computed
# as the cumulative sum of daily melt from October 1st to August 1st.
# Third, the cumulative sum of daily melt was normalized by the total
# annual melt (on August 1st ) to estimate the (daily) fraction of cumulative
# annual melt, which was then sampled on the dates of maximum SWE (FMmax)
# and April 1st (FMApr1). The August 1st end date was chosen to avoid rare
# cases where early-season (i.e., September) snowfall could impact estimates
# of the total annual melt and to ensure that any late-lying snow (almost
# always gone by August at sensor locations) was recorded as annual melt.

valid_water_yrs <- snow_data_snotel %>%
  filter(!is.na(snow_water_equivalent)) %>%
  arrange(date) %>%
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  group_by(water_year) %>%
  summarize(
    n_days = n(),
    .groups = "drop"
  ) %>%
  filter(n_days >= 365)

frac_pre_swe <- snow_data_snotel %>%
  arrange(date) %>%
  # one missing value but it's july and 0 long before and after
  mutate(
    melt = snow_water_equivalent - lag(snow_water_equivalent, 1),
    melt = ifelse(melt > 0 & !is.na(melt), 0, melt),
    year = lubridate::year(date),
    month = lubridate::month(date),
    doy = lubridate::yday(date),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  inner_join(., snow_phenology_snotel) %>%
  filter(water_year %in% valid_water_yrs$water_year) %>%
  mutate(bef = ifelse(doy < max_swe_doy | month >= 10, "bef", "after")) %>%
  # per keith's methods 0 out any melt in Aug, Sept
  # this doesn't actually matter here...
  mutate(melt = ifelse(month %in% c(8, 9), 0, melt)) %>%
  group_by(bef, water_year) %>%
  summarize(melt_acc = -1 * sum(melt), .groups = "drop") %>%
  pivot_wider(names_from = bef, values_from = melt_acc) %>%
  mutate(
    tot = after + bef,
    frac = bef / tot
  )


# calc fraction melted before april 1 instead?
prior_melt_apr_1 <- snow_data_snotel %>%
  arrange(date) %>%
  # one missing value but it's july and 0 long before and after
  mutate(
    snow_water_equivalent = ifelse(date == "2022-07-16", 0, snow_water_equivalent),
    melt = snow_water_equivalent - lag(snow_water_equivalent, 1),
    melt = ifelse(melt > 0 & !is.na(melt), 0, melt),
    year = lubridate::year(date),
    month = lubridate::month(date),
    doy = lubridate::yday(date),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  inner_join(., snow_phenology_snotel) %>%
  filter(water_year %in% valid_water_yrs$water_year) %>%
  mutate(bef = ifelse(month <= 3 | month >= 10, "bef", "after")) %>%
  # per keith's methods 0 out any melt in Aug, Sept
  # this doesn't actually matter here...
  mutate(melt = ifelse(month %in% c(8, 9), 0, melt)) %>%
  group_by(bef, water_year) %>%
  summarize(melt_acc = -1 * sum(melt)) %>%
  pivot_wider(names_from = bef, values_from = melt_acc) %>%
  mutate(
    tot = after + bef,
    frac = bef / tot
  )

#

# Calculate trends by site and season-------------------------------------------
trend_results_snowmelt_doy <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "last_snow_melt_doy", year_col = "year")
trend_results_max_swe_doy <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "max_swe_doy", year_col = "year")
trend_results_max_swe <- snow_phenology_snotel %>%
  perform_trend_analysis(., value_col = "max_swe", year_col = "year")
trend_results_max_frac_pre <- frac_pre_swe %>%
  perform_trend_analysis(., value_col = "frac", year_col = "water_year")
trend_results_prior_melt_apr_1 <- prior_melt_apr_1 %>%
  perform_trend_analysis(., value_col = "frac", year_col = "water_year")


# Create individual styled plots matching c1_d1_temps_annual format ---------------

# Helper function to create styled individual plots
create_individual_snow_plot <- function(data, trend_result, y_col, y_label, plot_title,
                                        slope_units, y_min = NULL, y_max = NULL) {
  # Prepare data
  plot_data <- data %>%
    mutate(
      season = "annual",
      site = "c1",
      year = year
    ) %>%
    left_join(trend_result %>% mutate(site = "c1", season = "annual"),
      by = c("site", "season")
    ) %>%
    mutate(
      significant = ifelse(mk_pvalue < 0.05, "Significant", "Non-significant"),
      significant = factor(significant, levels = c("Significant", "Non-significant"))
    )

  # Calculate trend line
  trend_line <- plot_data %>%
    filter(!is.na(ts_slope) & !is.na(ts_intercept)) %>%
    summarize(
      min_year = min(year, na.rm = TRUE),
      max_year = max(year, na.rm = TRUE),
      ts_slope = first(ts_slope),
      ts_intercept = first(ts_intercept),
      mk_pvalue = first(mk_pvalue),
      significant = first(significant),
      site = first(site)
    ) %>%
    mutate(
      y_start = ts_intercept + ts_slope * min_year,
      y_end = ts_intercept + ts_slope * max_year
    )

  # Create legend label with slope
  slope_value <- round(trend_line$ts_slope, 3)
  p_round <- round(as.numeric(trend_line$mk_pvalue), 3)
  legend_label <- paste0(ifelse(slope_value > 0, "+", ""), slope_value, " ", slope_units, " ; p-value: ", p_round)

  # Main plot
  main_plot <- ggplot(plot_data, aes(x = year, y = !!sym(y_col), color = site)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_line(linewidth = 1) +
    geom_segment(
      data = trend_line,
      aes(
        x = min_year, y = y_start,
        xend = max_year, yend = y_end,
        color = site,
        linetype = significant
      ),
      linewidth = 1.2
    ) +
    scale_color_manual(
      values = c("c1" = "#D55E00"),
      labels = c("c1" = legend_label),
      name = "Site"
    ) +
    scale_linetype_manual(
      values = c("Significant" = "solid", "Non-significant" = "dashed"),
      name = "Mann-Kendall\nSignificance"
    ) +
    scale_x_continuous(breaks = seq(1980, 2030, by = 10)) +
    labs(
      title = plot_title,
      x = "Year",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(
        order = 1, nrow = 1, override.aes = list(linewidth = 1.2), keywidth = 1,
        title.theme = element_text(size = 8),
        label.theme = element_text(size = 8)
      ),
      linetype = guide_legend(
        order = 2, nrow = 1, override.aes = list(linewidth = 1.2), keywidth = 2.5,
        title.theme = element_text(size = 8),
        label.theme = element_text(size = 8)
      )
    )

  # Extract site legend for positioning inside plot
  site_legend_plot <- ggplot(plot_data, aes(x = year, y = !!sym(y_col), color = site)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_line() +
    scale_color_manual(
      values = c("c1" = "#D55E00"),
      labels = c("c1" = legend_label),
      name = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.margin = margin(6, 6, 6, 6)
    ) +
    guides(color = guide_legend(
      override.aes = list(linewidth = 1.2, alpha = 0.6, size = 2.5),
      label.theme = element_text(size = 8)
    ))

  site_legend <- get_legend(site_legend_plot)

  # Extract significance legend for bottom
  sig_legend_plot <- ggplot(
    data.frame(x = 1:2, y = 1:2, sig = factor(c("Significant", "Non-significant"),
      levels = c("Significant", "Non-significant")
    )),
    aes(x = x, y = y, linetype = sig)
  ) +
    geom_line() +
    scale_linetype_manual(
      values = c("Significant" = "solid", "Non-significant" = "dashed"),
      name = "Mann-Kendall\nSignificance"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(linetype = guide_legend(
      nrow = 1, override.aes = list(linewidth = 1.2),
      keywidth = 2.5,
      title.theme = element_text(size = 8),
      label.theme = element_text(size = 8)
    ))

  sig_legend <- get_legend(sig_legend_plot)

  # Create plot without legend
  plot_no_legend <- main_plot + theme(legend.position = "none")

  # Add site legend inside plot
  # Auto-calculate positioning if not provided
  if (is.null(y_min)) y_min <- min(plot_data[[y_col]], na.rm = TRUE)
  if (is.null(y_max)) y_max <- max(plot_data[[y_col]], na.rm = TRUE)

  y_range <- y_max - y_min
  legend_ymin <- y_min + (y_range * 0.003)
  legend_ymax <- y_min + (y_range * 0.09)

  plot_with_legend <- plot_no_legend +
    annotation_custom(
      grob = site_legend,
      xmin = 2010, xmax = 2025,
      ymin = legend_ymin, ymax = legend_ymax
    )

  # Return both the plot with internal legend and the sig legend
  return(list(plot = plot_with_legend, sig_legend = sig_legend))
}

# Extract legend from one plot
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Create all individual plots
max_swe_doy_result <- create_individual_snow_plot(
  data = snow_phenology_snotel,
  trend_result = trend_results_max_swe_doy,
  y_col = "max_swe_doy",
  y_label = "Maximum SWE (DOY)",
  plot_title = "Timing of Maximum SWE",
  slope_units = "days per year"
)

snowmelt_doy_result <- create_individual_snow_plot(
  data = snow_phenology_snotel,
  trend_result = trend_results_snowmelt_doy,
  y_col = "last_snow_melt_doy",
  y_label = "Last snowmelt (DOY)",
  plot_title = "Timing of Last Snowmelt",
  slope_units = "days per year"
)

max_swe_result <- create_individual_snow_plot(
  data = snow_phenology_snotel,
  trend_result = trend_results_max_swe,
  y_col = "max_swe",
  y_label = "Maximum SWE (mm)",
  plot_title = "Maximum SWE",
  slope_units = "mm per year"
)

frac_pre_result <- create_individual_snow_plot(
  data = frac_pre_swe %>% rename(year = water_year) %>% ungroup(),
  trend_result = trend_results_max_frac_pre,
  y_col = "frac",
  y_label = "Fraction of melt\nbefore maximum SWE",
  plot_title = "Fraction of Winter Melt",
  slope_units = "% per year"
)


# Save all individual plots to figures folder
if (!is.null(dev.list())) dev.off() # Close any open devices

ggsave(
  filename = file.path(figures_dir, "max_swe_doy.jpg"),
  plot = grid.arrange(max_swe_doy_result$plot, max_swe_doy_result$sig_legend,
    nrow = 2, heights = c(5, 0.5)
  ),
  width = 16, height = 12, units = "in", dpi = 300,
  scale = 0.5
)

ggsave(
  filename = file.path(figures_dir, "last_snowmelt_doy.jpg"),
  plot = grid.arrange(snowmelt_doy_result$plot, snowmelt_doy_result$sig_legend,
    nrow = 2, heights = c(5, 0.5)
  ),
  width = 16, height = 12, units = "in", dpi = 300,
  scale = 0.5
)

ggsave(
  filename = file.path(figures_dir, "max_swe_mm.jpg"),
  plot = grid.arrange(max_swe_result$plot, max_swe_result$sig_legend,
    nrow = 2, heights = c(5, 0.5)
  ),
  width = 16, height = 12, units = "in", dpi = 300,
  scale = 0.5
)

ggsave(
  filename = file.path(figures_dir, "frac_melt_pre_max_swe.jpg"),
  plot = grid.arrange(frac_pre_result$plot, frac_pre_result$sig_legend,
    nrow = 2, heights = c(5, 0.5)
  ),
  width = 16, height = 12, units = "in", dpi = 300,
  scale = 0.5
)

print("All individual plots saved to figures folder!")

# Combined plot ----------------------------------------------------------------
if (!is.null(dev.list())) dev.off() # closes any open devices

# Create a temporary plot with legend just for extraction
temp_plot_with_legend <- ggplot(
  snow_phenology_snotel %>%
    mutate(season = "annual", site = "c1"),
  aes(x = year, y = max_swe, color = site)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line() +
  scale_color_manual(
    values = c("c1" = "#D55E00"),
    name = NULL,
    labels = c("c1" = "Subalpine")
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linewidth = 1,
      shape = 16,
      size = 1.5,
      alpha = 0.6
    )
  ))

shared_legend <- get_legend(temp_plot_with_legend)

# Create a temporary plot to extract significance legend
sig_legend_combined <- ggplot(
  data.frame(x = 1:3, y = 1:3, sig = factor(c("Significant", "Non-significant", "Insufficient data"),
    levels = c("Significant", "Non-significant", "Insufficient data")
  )),
  aes(x = x, y = y, linetype = sig)
) +
  geom_line() +
  scale_linetype_manual(
    values = c("Significant" = "solid", "Non-significant" = "dashed", "Insufficient data" = "dotted"),
    name = "Mann-Kendall Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(linetype = guide_legend(
    nrow = 3, override.aes = list(linewidth = 1.2),
    keywidth = 2.5,
    title.theme = element_text(size = 14),
    label.theme = element_text(size = 13)
  ))

sig_legend <- get_legend(sig_legend_combined)

# Combine both legends horizontally
combined_legends <- arrangeGrob(shared_legend, sig_legend,
  ncol = 2, widths = c(1, 1)
)

jpeg(file.path(figures_dir, "combined_plot_snow.jpg"), width = 12, height = 10, units = "in", res = 300)

# Create plots WITHOUT tables for the combined figure
grid.arrange(
  create_seasonal_plot(
    seasonal_data = (snow_phenology_snotel %>%
      mutate(season = "annual", site = "c1")),
    trend_results = (trend_results_max_swe_doy %>%
      mutate(site = "c1", season = "annual")),
    target_season = "annual",
    y_col = "max_swe_doy", y_label = "Maximum SWE (DOY)",
    plot_title = "Timing of Maximum SWE", slope_units = "days per year",
    include_table = FALSE
  ),
  create_seasonal_plot(
    seasonal_data = (snow_phenology_snotel %>%
      mutate(season = "annual", site = "c1")),
    trend_results = (trend_results_snowmelt_doy %>%
      mutate(site = "c1", season = "annual")),
    target_season = "annual",
    y_col = "last_snow_melt_doy", y_label = "Last snowmelt (DOY)",
    plot_title = "Timing of Last Snowmelt", slope_units = "days per year",
    include_table = FALSE
  ),
  create_seasonal_plot(
    seasonal_data = (snow_phenology_snotel %>%
      mutate(season = "annual", site = "c1")),
    trend_results = (trend_results_max_swe %>%
      mutate(site = "c1", season = "annual")),
    target_season = "annual",
    y_col = "max_swe", y_label = "Maximum SWE (mm)",
    plot_title = "Maximum SWE", slope_units = "mm per year",
    include_table = FALSE
  ),
  create_seasonal_plot(
    seasonal_data = (frac_pre_swe %>%
      mutate(season = "annual", site = "c1") %>%
      rename(year = water_year) %>%
      ungroup()),
    trend_results = (trend_results_max_frac_pre %>%
      mutate(site = "c1", season = "annual")),
    target_season = "annual",
    y_col = "frac", y_label = "Fraction of melt\nbefore maximum SWE",
    plot_title = "Fraction of Winter Melt", slope_units = "% per year",
    include_table = FALSE
  ),
  combined_legends,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
  heights = c(1, 1, 0.15)
)
dev.off()

print("Combined plot saved to figures folder!")
