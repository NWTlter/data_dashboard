# data_dashboard
Niwot LTER data dashboard

## Repository Structure

```
data_dashboard/
├── README.md                           # This file - project overview and documentation
├── LICENSE                            # Project license
├── data_dashboard.Rproj               # RStudio project file
│
├── c1_d1_sdl_temp_ppt/               # Temperature and precipitation analysis for C1, D1, and SDL sites
│   ├── plot_c1_d1_sdl_temp_ppt_infilled.R  # Main script: analyzes temperature and precipitation trends 
│   │                                        # by season, performs Mann-Kendall trend tests
│   ├── data/                          # Climate data from multiple NWT LTER datasets
│   └── figures/                       # Generated climate trend visualizations
│       ├── c1_d1_summer_mean_temp.jpg # Summer temperature trends comparison
│       ├── combined_plot_ppt.jpg      # Combined precipitation analysis plot
│       ├── combined_plot_temp.jpg     # Combined temperature analysis plot
│       ├── d1_ppt_anom_by_season.png  # D1 precipitation anomalies by season
│       └── d1_temp_anom_by_season.png # D1 temperature anomalies by season
│
├── c1_snow/                          # Snow analysis for C1 site
│   ├── plot_frac_winter_melt_trends.R # Analyzes winter snow melt fraction trends using SNOTEL data
│   └── figures/
│       └── combined_plot_snow.jpg    # Snow trend analysis visualization
│
├── noaa_gml/                         # NOAA Global Monitoring Laboratory data analysis
│   ├── download_noaa_gml.R           # Downloads atmospheric composition data from NOAA GML
│   ├── plot_noaa_co2_consortium.R    # Creates CO2 trend plots from Niwot Ridge measurements
│   └── figures/
│       └── noaa_co2.png              # CO2 trend visualization for Niwot Ridge
│
├── sdl_moisture/                     # Soil moisture analysis for Saddle site
│   ├── plot_sdl_soil_moist.R        # Analyzes soil moisture and temperature patterns (NWT package 405)
│   └── figures/                      # Soil moisture analysis plots (generated)
│
├── sdl_snow/                         # Snow analysis for Saddle site
│   ├── est_sdl_snowmelt_doy.R        # Estimates snowmelt day-of-year using Bayesian interval 
│   │                                 # censored models from snow survey data
│   └── data_deriv/                   # Derived data products
│       └── snowmelt_est_adj.csv      # Estimated snowmelt dates with adjustments
│
└── sdl_veg/                          # Vegetation composition analysis for Saddle site
    └── figures/                         # Vegetation analysis plots
        ├── sdl_veg.jpg              # Growth form cover anomalies by year
```

## Project Overview

This repository contains data analysis scripts and visualizations for the Niwot LTER data dashboard. The project analyzes long-term environmental trends across multiple research sites in the Niwot Ridge ecosystem, including:

- **Climate trends**: Temperature and precipitation patterns at C1, D1, and Saddle (SDL) sites
- **Snow dynamics**: Winter snow accumulation, melt timing, and trends
- **Atmospheric composition**: CO2 measurements from NOAA Global Monitoring Laboratory
- **Soil conditions**: Soil moisture patterns
- **Vegetation**: Plant community composition changes by growth form

Each analysis folder contains:
- **R scripts**: Data processing, analysis, and visualization code
- **data/**: Raw and processed datasets, often downloaded from EDI (Environmental Data Initiative)
- **figures/**: Generated plots and visualizations for the dashboard

## Data Sources

- **EDI (Environmental Data Initiative)**: Primary repository for NWT LTER datasets
- **NOAA GML**: Atmospheric composition measurements from Niwot Ridge
- **SNOTEL**: Snow telemetry data for snow trend analysis
