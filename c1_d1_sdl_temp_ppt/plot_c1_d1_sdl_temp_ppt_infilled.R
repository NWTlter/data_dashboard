# plot c1/d1 ppt trends by quarter
# and/or by winter/summer
# SCE 11 July 2025

# -- SETUP ------
# clean environment, load needed libraries, modify default settings
rm(list = ls())
library(tidyverse)
library(lubridate)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c("NP", "NA", NA, "NaN", NaN, ".")

library(tidyverse)
library(EDIutils) # Handy tools for interacting with EDI's API
library(trend) # for Mann-Kendall and Theil-Sen tests
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(trend)

# only need to download once
download_data <- FALSE

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
  n_sites <- length(unique(plot_data$site))
  site_colors <- rainbow(n_sites)
  names(site_colors) <- unique(plot_data$site)

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
    )

  # Create the plot title
  if (is.null(plot_title)) {
    plot_title <- paste(stringr::str_to_title(target_season), "Trends by Site")
  }

  # Create main plot using tidy evaluation for dynamic y column
  p <- ggplot(plot_data, aes(x = year, y = !!sym(y_col), color = site)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_segment(
      data = trend_lines,
      aes(
        x = min_year, y = y_start,
        xend = max_year, yend = y_end,
        color = site, linetype = significant
      ),
      size = 1.2
    ) +
    scale_color_manual(values = site_colors, name = "Site") +
    scale_linetype_manual(
      values = c(
        "Significant" = "solid",
        "Non-significant" = "dashed",
        "Insufficient data" = "dotted"
      ),
      name = "Mann-Kendall\nSignificance"
    ) +
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

  # Create inset table with slope information
  slope_data <- plot_data %>%
    select(site, ts_slope, mk_pvalue, mk_trend) %>%
    distinct() %>%
    mutate(
      slope_round = round(ts_slope, 3),
      p_round = round(mk_pvalue, 3),
      significance = ifelse(mk_pvalue < 0.05, "*", ""),
      display_text = paste0(slope_round, significance)
    )

  # Create inset grob with dynamic slope units
  inset_data <- data.frame(
    Site = slope_data$site,
    stringsAsFactors = FALSE
  )
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
    color <- site_colors[site_name]
    table_grob <- gtable_add_grob(table_grob,
      rectGrob(gp = gpar(fill = color, alpha = 0.3, col = NA)),
      t = i + 1, l = 1, r = 1
    )
  }

  # Add the table as annotation - use dynamic y column for positioning
  y_values <- plot_data[[y_col]]
  p_with_inset <- p +
    annotation_custom(table_grob,
      xmin = max(plot_data$year) - (max(plot_data$year) - min(plot_data$year)) * 0.35,
      xmax = max(plot_data$year) - (max(plot_data$year) - min(plot_data$year)) * 0.05,
      ymin = max(y_values, na.rm = TRUE) - (max(y_values, na.rm = TRUE) - min(y_values, na.rm = TRUE)) * 0.4,
      ymax = max(y_values, na.rm = TRUE) - (max(y_values, na.rm = TRUE) - min(y_values, na.rm = TRUE)) * 0.05
    )

  return(p_with_inset)
}



# download data -----------------------------------------------------------
# note if you have already downloaded SOME data the read_data_package_archive
# function will bork as it doesn't want to overwrite, so clear your /data
# directories and then rerun
if (download_data) {
  # download the data from EDI
  # 184 is c1 ppt
  # 186 is d1 ppt
  # 315 is the sdl ppt
  # 314 is sdl temp
  # 185 is c1 temp
  # 187 is d1 temp

  scope <- "knb-lter-nwt" # Niwot scope

  # note the overwrite argument does not work so clear out any existing
  # copies before running this
  for (id in c("184", "186", "315",
    "314", "185", "187"
  )) {
    # ask EDI to tell me what the most current version is
    revision <- list_data_package_revisions(scope, id, filter = "newest")

    # display current version - > this is referred to as the "packageID"
    packageID <- paste(scope, id, revision, sep = ".")

    # download the data
    read_data_package_archive(packageID, path = "c1_d1_sdl_temp_ppt/data")
    print(read_data_package_citation(packageID))
  }

  # update the below so you remember to cite it correctly
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for C1 chart recorder, 1952 - ongoing. ver 8. Environmental Data Initiative. https://doi.org/10.6073/pasta/5f896d7d1eb19649ee2ff9fca7166160. Accessed 2025-07-09."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for D1 chart recorder, 1952 - ongoing. ver 6. Environmental Data Initiative. https://doi.org/10.6073/pasta/02b125f7d0dc87a6d40e6ca7001471ce. Accessed 2025-07-09."
  # [1] "White, C., J. Morse, K. Chowanski, T. Kittel, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily precipitation data for Saddle, 1981 - ongoing. ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/8caa444276521b797e4a16485dc0e137. Accessed 2025-07-11."
  # [1] "White, C., J. Morse, H. Brandes, K. Chowanski, T. Kittel, M. Losleben, and M. Moore. 2025. Homogenized, gap-filled, daily air temperature data for Saddle, 1986 - ongoing. ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/49a2171ec028c4ddef298752bc9ebe8d. Accessed 2025-07-11."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for C1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/0b4f8747ea72258e46681bb511c262f3. Accessed 2025-07-11."
  # [1] "Kittel, T., C. White, M. Hartman, K. Chowanski, T. Ackerman, M. Williams, M. Losleben, and M. Moore. 2025. Infilled daily air temperature data for D1 chart recorder, 1952 - ongoing. ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/b05689181c21ac40d35b6c1c01e2f8a5. Accessed 2025-07-11."
  # unzip all zip files in the directory
  # overwrites the manifests but don't really need them.
  for (fname in list.files("c1_d1_sdl_temp_ppt/data",
    pattern = "knb-lter.*zip",
    full.names = TRUE
  )) {
    unzip(zipfile = fname, exdir = "c1_d1_sdl_temp_ppt/data/")
  }
}

# read and munge ppt and temp --------------------------------------------------

ppt_seasonal <- read_csv("c1_d1_sdl_temp_ppt/data/sdl_daily_precip_gapfilled.cw.data.csv") %>%
  bind_rows(., read_csv("c1_d1_sdl_temp_ppt/data/d1_infilled_precip_daily.tk.data.csv")) %>%
  bind_rows(., read_csv("c1_d1_sdl_temp_ppt/data/c1_infilled_precip_daily.tk.data.csv")) %>%
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
  filter(!(site == "SDL" & year < 1990))

temp_seasonal <- read_csv("c1_d1_sdl_temp_ppt/data/d1_infilled_temp_daily.tk.data.csv",
  guess_max = 100000
) %>%
  bind_rows(., read_csv("c1_d1_sdl_temp_ppt/data/c1_infilled_temp_daily.tk.data.csv",
    guess_max = 100000
  )) %>%
  bind_rows(., read_csv("c1_d1_sdl_temp_ppt/data/sdl_daily_airtemp_gapfilled.cw.data.csv",
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

# Calculate trends by site and season-------------------------------------------
trend_results_ppt <- ppt_seasonal %>%
  group_by(site, season) %>%
  do(perform_trend_analysis(., value_col = "tot_ppt", year_col = "year")) %>%
  ungroup()

trend_results_temp <- temp_seasonal %>%
  group_by(site, season) %>%
  do(perform_trend_analysis(., value_col = "mean_temp", year_col = "year")) %>%
  ungroup()

# Calculate trends by site and season-------------------------------------------
# Generate individual plots
summer_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "summer")
winter_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "winter")

summer_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "summer",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = NULL, slope_units = "degC/yr"
)
spring_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "spring",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = NULL, slope_units = "degC/yr"
)

fall_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "fall",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = NULL, slope_units = "degC/yr"
)

winter_plot_temp <- create_seasonal_plot(temp_seasonal, trend_results_temp,
  target_season = "winter",
  y_col = "mean_temp", y_label = "Temp (deg C)",
  plot_title = NULL, slope_units = "degC/yr"
)

summer_plot_temp <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "summer")
winter_plot_ppt <- create_seasonal_plot(ppt_seasonal, trend_results_ppt, "winter")


# combined plots ---------------------------------------------------------------

# Extract legend from one plot
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

shared_legend <- get_legend(spring_plot_temp)


combined_plot_ppt <- grid.arrange(winter_plot_ppt + theme(legend.position = "none"),
  summer_plot_ppt + theme(legend.position = "none"),
  shared_legend,
  ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1, 2), c(5, 5)),
  heights = c(1, 0.1)
)


combined_plot_temp <- grid.arrange(winter_plot_temp + theme(legend.position = "none"),
  spring_plot_temp + theme(legend.position = "none"),
  summer_plot_temp + theme(legend.position = "none"),
  fall_plot_temp + theme(legend.position = "none"),
  shared_legend,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
  heights = c(1, 1, 0.1)
)


