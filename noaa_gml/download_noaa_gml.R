# NOAA Data Download Automation Script
# This script automates downloading data files from the NOAA GML data page
# It saves a local copy of all the NOAA GML data in case of future changes 
# to the NOAA ftp

# Required libraries
library(httr)
library(dplyr)
library(stringr)
library (rvest)

# funtion to download using direct HTTP requests
download_noaa_data_http <- function(url = "https://gml.noaa.gov/data/data.php?site=NWR&perpage=200",
                                    download_dir = "noaa_gml") {
  
  # Create download directory
  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }
  
  # Read the page HTML
  cat("Fetching page content...\n")
  page <- read_html(url)
  
  # Extract all download links
  download_links <- page %>%
    html_nodes("a[href]") %>%
    html_attr("href") %>%
    .[str_detect(., "\\.(zip|txt|dat|csv)$|download|dataset.php")]
  
  download_links <- ifelse(str_detect(download_links, "^http"), 
                           download_links,
                           paste0("https://gml.noaa.gov/data/", download_links))
  
  cat(sprintf("Found %d download links\n", length(download_links)))
  
  # Download each file
  downloaded_files <- character()
  
  for (i in seq_along(download_links)) {
    tryCatch({
      download_url <- download_links[i]
      #https://gml.noaa.gov/data/data.php?site=NWR&perpage=200
      #https://gml.noaa.gov/data/dataset.php?item=nwr-CCl4-cats-daily # this is what you want
      
      download_file <- read_html(download_url)%>%
        html_nodes("a[href]") %>%
        html_attr("href") %>%
        .[str_detect(., "\\.(zip|txt|dat|csv)$|download|dataset.php")]
  
        filename <- basename(download_file)
      
      local_path <- file.path(download_dir, filename)
      
      cat(sprintf("Downloading %d/%d: %s\n", i, length(download_links), filename))
      
      download_url <- paste0('https://gml.noaa.gov', download_file)
      response <- GET(download_url)

      if (status_code(response) == 200) {
        writeBin(content(response, "raw"), local_path)
        downloaded_files <- c(downloaded_files, local_path)
        cat(sprintf("Successfully downloaded: %s\n", filename))
      } else {
        cat(sprintf("Failed to download: %s (Status: %d)\n", download_url, status_code(response)))
      }
      
      Sys.sleep(1)  # Be respectful to the server
      
    }, error = function(e) {
      cat(sprintf("Error downloading file %d: %s\n", i, e$message))
    })
  }
  
  cat(sprintf("\nHTTP download completed. %d files downloaded to: %s\n", 
              length(downloaded_files), download_dir))
  
  return(downloaded_files)
}

download_noaa_data_http(url = "https://gml.noaa.gov/data/data.php?site=NWR&perpage=200",
                                    download_dir = "noaa_gml/data")


